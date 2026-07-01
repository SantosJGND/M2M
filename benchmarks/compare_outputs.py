#!/usr/bin/env python3
"""Compare Rust vs Python output TSV files with approximate numeric tolerance."""

import sys
import os


def read_tsv(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    data = [l.split("\t") for l in lines[1:]]
    return header, data


def compare_tsvs(rust_path, python_path, atol=1e-6):
    if not os.path.exists(rust_path):
        return f"MISSING: {rust_path}"
    if not os.path.exists(python_path):
        return f"MISSING: {python_path}"

    rust_header, rust_data = read_tsv(rust_path)
    py_header, py_data = read_tsv(python_path)

    if rust_header != py_header:
        return f"HEADER MISMATCH: {rust_header} vs {py_header}"

    if len(rust_data) != len(py_data):
        return f"ROW COUNT MISMATCH: {len(rust_data)} vs {len(py_data)}"

    n_rows = len(rust_data)
    n_cols = len(rust_header)

    max_diff = 0.0
    mismatch_count = 0
    for i in range(n_rows):
        if len(rust_data[i]) != len(py_data[i]):
            return f"COL COUNT MISMATCH row {i}: {len(rust_data[i])} vs {len(py_data[i])}"
        for j in range(n_cols):
            try:
                rv = float(rust_data[i][j])
                pv = float(py_data[i][j])
            except ValueError:
                if rust_data[i][j] != py_data[i][j]:
                    mismatch_count += 1
                continue
            diff = abs(rv - pv)
            if diff > max_diff:
                max_diff = diff
            if diff > atol:
                mismatch_count += 1

    if mismatch_count == 0:
        return f"PASS (max_diff={max_diff:.2e})"
    else:
        return f"FAIL: {mismatch_count} mismatches (max_diff={max_diff:.2e}, tol={atol})"


def main():
    rust_dir = sys.argv[1] if len(sys.argv) > 1 else "test_output"
    python_dir = sys.argv[2] if len(sys.argv) > 2 else "python_output"

    files = ["distance_matrix.tsv", "presence_absence_matrix.tsv"]

    all_pass = True
    for fname in files:
        rust_path = os.path.join(rust_dir, fname)
        py_path = os.path.join(python_dir, fname)
        result = compare_tsvs(rust_path, py_path)
        print(f"{fname}: {result}")
        if result.startswith("FAIL") or result.startswith("MISSING"):
            all_pass = False

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
