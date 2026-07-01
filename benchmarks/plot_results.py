#!/usr/bin/env python3
"""Generate benchmark comparison plot from Rust and Python results."""

import json
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_results(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def main():
    rust_path = "benchmarks/rust_results.json"
    python_path = "benchmarks/python_results.json"

    rust_data = load_results(rust_path)
    python_data = load_results(python_path)

    if rust_data is None and python_data is None:
        print("ERROR: No benchmark results found.")
        print("Run 'cargo bench' (Rust) and 'python benchmarks/run_bench.py' (Python) first.")
        sys.exit(1)

    # Build combined results
    combined = {"metadata": {}}
    stages_order = ["parse_bam_files", "to_matrix_presence_absence",
                    "to_shared_mut_matrix", "standardize_by_files_total",
                    "full_pipeline"]

    if rust_data:
        combined["metadata"].update(rust_data.get("metadata", {}))
    if python_data:
        combined["metadata"].update(python_data.get("metadata", {}))

    combined["rust"] = rust_data.get("results", {}) if rust_data else {}
    combined["python"] = python_data.get("results", {}) if python_data else {}

    # Filter to stages present in both
    stages = [s for s in stages_order
              if s in combined.get("rust", {}) or s in combined.get("python", {})]

    if not stages:
        print("ERROR: No benchmark stages to plot.")
        sys.exit(1)

    # Build plot
    fig, ax = plt.subplots(figsize=(12, 7))

    x = np.arange(len(stages))
    w = 0.35

    rust_means = []
    rust_stds = []
    python_means = []
    python_stds = []

    for s in stages:
        r = combined.get("rust", {}).get(s, {})
        p = combined.get("python", {}).get(s, {})
        rust_means.append(r.get("mean_ms", 0))
        rust_stds.append(r.get("std_ms", 0))
        python_means.append(p.get("mean_ms", 0))
        python_stds.append(p.get("std_ms", 0))

    rust_bars = ax.bar(x - w/2, rust_means, w, yerr=rust_stds,
                       label="Rust", color="#e76f51", capsize=4, ec="black", lw=0.5)
    python_bars = ax.bar(x + w/2, python_means, w, yerr=python_stds,
                         label="Python", color="#264653", capsize=4, ec="black", lw=0.5)

    # Annotate speedup ratios
    for i, (rm, pm) in enumerate(zip(rust_means, python_means)):
        if rm > 0 and pm > 0:
            if pm >= rm:
                ratio = pm / rm
                label = f"{ratio:.1f}×"
                y_pos = max(rm, pm) + max(rust_stds[i], python_stds[i]) + (max(rust_means + python_means) * 0.02)
                ax.annotate(label, (x[i], y_pos), ha="center", va="bottom",
                            fontsize=9, fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.7))

    ax.set_xticks(x)
    ax.set_xticklabels(stages, rotation=20, ha="right", fontsize=10)
    ax.set_ylabel("Time (ms)", fontsize=12)
    ax.set_xlabel("Benchmark Stage", fontsize=12)

    meta = combined.get("metadata", {})
    title = f"Rust vs Python Benchmark\n"
    title += f"Data: {meta.get('test_data', '?')}  |  "
    title += f"Files: {meta.get('num_files', '?')}  |  "
    title += f"max_reads: {meta.get('max_reads', '?')}  |  "
    title += f"freq_threshold: {meta.get('frequency_threshold', '?')}"
    ax.set_title(title, fontsize=13, fontweight="bold")

    ax.legend(fontsize=11)
    ax.grid(axis="y", alpha=0.3)

    # Use log scale if large disparity
    all_vals = rust_means + python_means
    if max(all_vals) / max(min(all_vals), 1) > 20:
        ax.set_yscale("log")
        ax.set_ylabel("Time (ms, log scale)", fontsize=12)

    plt.tight_layout()

    os.makedirs("benchmarks", exist_ok=True)
    out_path = "benchmarks/benchmark_comparison.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Plot saved to {out_path}")

    # Write summary table
    table_path = "benchmarks/results_table.tsv"
    with open(table_path, "w") as f:
        f.write("Stage\tRust (ms)\tPython (ms)\tSpeedup\n")
        for i, s in enumerate(stages):
            rm = rust_means[i]
            pm = python_means[i]
            rs = rust_stds[i]
            ps = python_stds[i]
            if rm > 0 and pm > 0:
                ratio = f"{pm/rm:.1f}×" if pm >= rm else f"{rm/pm:.1f}× (Python faster)"
            else:
                ratio = "N/A"
            rust_str = f"{rm:.1f} ± {rs:.1f}" if rm > 0 else "N/A"
            py_str = f"{pm:.1f} ± {ps:.1f}" if pm > 0 else "N/A"
            f.write(f"{s}\t{rust_str}\t{py_str}\t{ratio}\n")
    print(f"Summary table saved to {table_path}")

    # Generate combined results JSON
    combined_path = "benchmarks/results.json"
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"Combined results saved to {combined_path}")


if __name__ == "__main__":
    main()
