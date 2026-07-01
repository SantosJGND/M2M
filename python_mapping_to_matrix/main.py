import os
import sys
import time

from .bam_parser import MultipleSamParser
from .matrix_ops import (
    to_shared_mut_matrix,
    to_matrix_presence_absence,
    standardize_by_files_total,
    invert_matrix,
    transpose_matrix,
)
from .cli import parse_args


def find_bam_files(directory):
    return [
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".bam")
    ]


def count_mapped_reads(files):
    import pysam
    counts = []
    for f in files:
        count = 0
        sam = pysam.AlignmentFile(f, "rb")
        for read in sam:
            if read.flag != 4:
                count += 1
        sam.close()
        counts.append(count)
    return counts


def write_distance_matrix(matrix, files, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "distance_matrix.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(files) + "\n")
        for i, row in enumerate(matrix):
            fh.write(files[i] + "\t" + "\t".join(f"{v:.6f}" for v in row) + "\n")


def write_presence_absence_matrix(matrix, files, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "presence_absence_matrix.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(files) + "\n")
        for row in matrix:
            fh.write("\t".join(str(v) for v in row) + "\n")


def main():
    args = parse_args()
    bam_files = []

    if args.files and args.input_directory:
        print("Please provide either a list of BAM files or an input directory, not both.")
        sys.exit(1)
    if not args.files and not args.input_directory:
        print("Please provide a list of BAM files or an input directory.")
        sys.exit(1)

    if args.files:
        bam_files = [f.strip() for f in args.files.split(",") if f.strip()]
    else:
        bam_files = find_bam_files(args.input_directory)
        if not bam_files:
            print("No bam files found in directory. Exiting.")
            sys.exit(1)

    bam_files = [f for f in bam_files if os.path.exists(f) and f.endswith(".bam")]
    if not bam_files:
        print("No valid BAM files found. Exiting.")
        sys.exit(1)

    parser = MultipleSamParser(
        frequency_threshold=args.frequency_threshold,
        max_reads=args.max_reads,
    )

    file_counts = count_mapped_reads(bam_files)
    print(f"[INFO] Mapped read counts per file: {file_counts}")

    parser.parse_files(bam_files, file_counts)
    print("[INFO] BAM files parsed (Phase 1).")
    parser.phase1_filter_reads_by_frequency()
    print("[INFO] Reads filtered by frequency.")
    parser.sample_reads()
    print("[INFO] Reads sampled (limiting to max_reads for Phase 2).")
    parser.guard_file_coverage()
    print("[INFO] File coverage guarded.")
    parser.phase2_parse_files()
    print("[INFO] BAM files parsed (Phase 2 - sampled reads only).")
    parser.filter_samples_to_keep(False)
    print("[INFO] Samples filtered.")

    if not parser.reads_to_keep:
        print("No reads to process after filtering. Exiting.")
        sys.exit(1)

    matrix = to_shared_mut_matrix(parser)
    print(f"[INFO] Distance matrix built with {len(matrix)} samples.")

    if len(matrix) == 0 or len(matrix[0]) == 0:
        print("Distance matrix is empty. Exiting.")
        sys.exit(1)

    shared_matrix = standardize_by_files_total(parser, matrix)
    inverted_distance_matrix = invert_matrix(shared_matrix)

    write_distance_matrix(inverted_distance_matrix, parser.files, args.output_directory)

    pa_matrix = to_matrix_presence_absence(parser)
    print(f"[INFO] Presence/absence matrix built with {pa_matrix.shape[0]} reads and {pa_matrix.shape[1]} samples.")
    tpa = transpose_matrix(pa_matrix)
    write_presence_absence_matrix(tpa, parser.files, args.output_directory)

    print("[INFO] Reports written to output directory.")


if __name__ == "__main__":
    main()
