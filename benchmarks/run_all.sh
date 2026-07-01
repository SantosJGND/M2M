#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

echo "=== mapping_to_matrix Benchmark Suite ==="
echo "Date: $(date)"
echo "Repo: $REPO_DIR"
echo ""

cd "$REPO_DIR"

# Ensure test data exists
if [ ! -d "test_data_larger" ] || [ -z "$(ls -A test_data_larger/*.bam 2>/dev/null)" ]; then
    echo "ERROR: test_data_larger directory missing or empty"
    exit 1
fi

# 1. Run Rust benchmarks
echo "----------------------------------------"
echo "  Running Rust benchmarks (cargo bench)..."
echo "----------------------------------------"
if command -v cargo &>/dev/null; then
    cargo bench 2>&1 | tee benchmarks/rust_bench_output.txt || echo "  WARN: cargo bench failed"
else
    echo "  WARN: cargo not found, skipping Rust benchmarks"
fi

# 2. Run Python benchmarks
echo ""
echo "----------------------------------------"
echo "  Running Python benchmarks..."
echo "----------------------------------------"
python3 "$SCRIPT_DIR/run_bench.py"

# 3. Parse Rust Criterion output into JSON (if available)
if [ -f "benchmarks/rust_bench_output.txt" ]; then
    echo ""
    echo "  Parsing Rust benchmark results..."
    python3 -c "
import json, re

results = {}
with open('benchmarks/rust_bench_output.txt') as f:
    text = f.read()

# Normalize multi-line Criterion output (name on separate line from time:)
lines = text.split(chr(10))
normalized = []
for line in lines:
    stripped = line.strip()
    if stripped.startswith('time:'):
        if normalized and not normalized[-1].strip().startswith('time:'):
            normalized[-1] += ' ' + stripped
        else:
            normalized.append(stripped)
    else:
        normalized.append(line)
full_text = chr(10).join(normalized)

pattern = r'(\w[\w_]+)\s+time:\s+\[([\d.]+)\s+ms\s+([\d.]+)\s+ms\s+([\d.]+)\s+ms\]'
for m in re.finditer(pattern, full_text):
    name = m.group(1)
    lo = float(m.group(2))
    mid = float(m.group(3))
    hi = float(m.group(4))
    std = (hi - lo) / 2.0
    results[name] = {'mean_ms': mid, 'std_ms': std, 'min_ms': lo, 'max_ms': hi}

if results:
    data = {
        'metadata': {
            'date': '2026-07-01',
            'test_data': 'test_data_larger',
            'language': 'rust',
        },
        'results': results,
    }
    with open('benchmarks/rust_results.json', 'w') as f:
        json.dump(data, f, indent=2)
    print('  Rust results written to benchmarks/rust_results.json')
    for name, r in results.items():
        print(f'    {name}: {r[\"mean_ms\"]:.3f} ms')
"
fi

# 4. Generate plot and summary table
echo ""
echo "----------------------------------------"
echo "  Generating plot and summary..."
echo "----------------------------------------"
python3 "$SCRIPT_DIR/plot_results.py"

echo ""
echo "=== Benchmark complete ==="
echo "Output files:"
echo "  benchmarks/benchmark_comparison.png  (plot)"
echo "  benchmarks/results_table.tsv         (table)"
echo "  benchmarks/results.json              (combined results)"
