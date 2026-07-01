#!/usr/bin/env python3
"""Per-stage Python benchmarks matching Rust Criterion bench structure.
Each bench times the ENTIRE operation (parse+filter+sample+stage) from scratch,
exactly matching how Criterion wraps everything in b.iter()."""

import sys
import os
import time
import json
import statistics

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from python_mapping_to_matrix.bam_parser import MultipleSamParser
from python_mapping_to_matrix.matrix_ops import (
    to_shared_mut_matrix,
    to_matrix_presence_absence,
    standardize_by_files_total,
    invert_matrix,
)
from python_mapping_to_matrix.main import find_bam_files

TEST_DATA_DIR = "test_data_larger"
FREQ_THRESHOLD = 0.1
MAX_READS = 500000
NUM_ITERATIONS = 30


def bench_parse_bam_files():
    times = []
    for _ in range(NUM_ITERATIONS):
        t0 = time.perf_counter()
        bam_files = find_bam_files(TEST_DATA_DIR)
        parser = MultipleSamParser(FREQ_THRESHOLD, MAX_READS)
        parser.parse_bam_files(bam_files)
        parser.reset_reads()
        parser.filter_reads_to_keep()
        parser.sample_reads()
        parser.filter_samples_to_keep(False)
        _ = (parser.read_count, len(parser.files))
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000)
    return {
        "mean_ms": statistics.mean(times),
        "std_ms": statistics.stdev(times) if len(times) > 1 else 0.0,
    }


def bench_to_matrix_presence_absence():
    times = []
    for _ in range(NUM_ITERATIONS):
        t0 = time.perf_counter()
        bam_files = find_bam_files(TEST_DATA_DIR)
        parser = MultipleSamParser(FREQ_THRESHOLD, MAX_READS)
        parser.parse_bam_files(bam_files)
        parser.reset_reads()
        parser.filter_reads_to_keep()
        parser.sample_reads()
        parser.filter_samples_to_keep(False)
        _ = to_matrix_presence_absence(parser)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000)
    return {
        "mean_ms": statistics.mean(times),
        "std_ms": statistics.stdev(times) if len(times) > 1 else 0.0,
    }


def bench_to_shared_mut_matrix():
    times = []
    for _ in range(NUM_ITERATIONS):
        t0 = time.perf_counter()
        bam_files = find_bam_files(TEST_DATA_DIR)
        parser = MultipleSamParser(FREQ_THRESHOLD, MAX_READS)
        parser.parse_bam_files(bam_files)
        parser.reset_reads()
        parser.filter_reads_to_keep()
        parser.sample_reads()
        parser.filter_samples_to_keep(False)
        _ = to_shared_mut_matrix(parser)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000)
    return {
        "mean_ms": statistics.mean(times),
        "std_ms": statistics.stdev(times) if len(times) > 1 else 0.0,
    }


def bench_standardize_by_files_total():
    times = []
    for _ in range(NUM_ITERATIONS):
        t0 = time.perf_counter()
        bam_files = find_bam_files(TEST_DATA_DIR)
        parser = MultipleSamParser(FREQ_THRESHOLD, MAX_READS)
        parser.parse_bam_files(bam_files)
        parser.reset_reads()
        parser.filter_reads_to_keep()
        parser.sample_reads()
        parser.filter_samples_to_keep(False)
        matrix = to_shared_mut_matrix(parser)
        _ = standardize_by_files_total(parser, matrix)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000)
    return {
        "mean_ms": statistics.mean(times),
        "std_ms": statistics.stdev(times) if len(times) > 1 else 0.0,
    }


def bench_full_pipeline():
    times = []
    for _ in range(NUM_ITERATIONS):
        t0 = time.perf_counter()
        bam_files = find_bam_files(TEST_DATA_DIR)
        parser = MultipleSamParser(FREQ_THRESHOLD, MAX_READS)
        parser.parse_bam_files(bam_files)
        parser.reset_reads()
        parser.filter_reads_to_keep()
        parser.sample_reads()
        parser.filter_samples_to_keep(False)
        if parser.reads_to_keep:
            matrix = to_shared_mut_matrix(parser)
            shared_matrix = standardize_by_files_total(parser, matrix)
            _ = invert_matrix(shared_matrix)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000)
    return {
        "mean_ms": statistics.mean(times),
        "std_ms": statistics.stdev(times) if len(times) > 1 else 0.0,
    }


def main():
    print(f"Python benchmarks ({NUM_ITERATIONS} iter each, full setup in timed section)")
    print(f"Data: {TEST_DATA_DIR}  max_reads={MAX_READS}  freq={FREQ_THRESHOLD}")
    print()

    results = {}

    for name, fn in [
        ("parse_bam_files", bench_parse_bam_files),
        ("to_matrix_presence_absence", bench_to_matrix_presence_absence),
        ("to_shared_mut_matrix", bench_to_shared_mut_matrix),
        ("standardize_by_files_total", bench_standardize_by_files_total),
        ("full_pipeline", bench_full_pipeline),
    ]:
        print(f"  {name}...", end=" ", flush=True)
        results[name] = fn()
        print(f"{results[name]['mean_ms']:.3f} ms")

    print()

    out = {
        "metadata": {
            "date": "2026-07-01",
            "test_data": TEST_DATA_DIR,
            "num_files": len(find_bam_files(TEST_DATA_DIR)),
            "max_reads": MAX_READS,
            "frequency_threshold": FREQ_THRESHOLD,
            "language": "python",
            "iterations": NUM_ITERATIONS,
        },
        "results": results,
    }

    os.makedirs("benchmarks", exist_ok=True)
    with open("benchmarks/python_results.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Results written to benchmarks/python_results.json")


if __name__ == "__main__":
    main()
