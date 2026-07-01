import pysam
import random
from collections import defaultdict


class RecordSummary:
    def __init__(self, nmatches: int, length: int):
        self.nmatches = nmatches
        self.length = length

    def percent_identity(self) -> float:
        return self.nmatches / self.length if self.length > 0 else 0.0


def matches_from_cigar(cigar_tuples) -> int:
    total = 0
    for op, length in cigar_tuples:
        if op == 0:
            total += length
    return total


class MultipleSamParser:
    def __init__(self, frequency_threshold: float = 0.1, max_reads: int = 500000):
        self.read_hash = {}
        self.read_names = []
        self.read_names_set = set()
        self.files = []
        self.reads_to_keep = []
        self.read_count = 0
        self.frequency_threshold = frequency_threshold
        self.max_reads = max_reads
        self.read_presence = {}
        self.phase1_file_counts = {}

    def parse_bam_files(self, bam_files):
        for bam_file in bam_files:
            file_reads = {}
            read_names_in_file = set()
            try:
                sam = pysam.AlignmentFile(bam_file, "rb")
            except Exception as e:
                print(f"[WARN] Failed to open {bam_file}: {e}")
                continue

            for read in sam:
                if read.flag == 4:
                    continue
                read_name = read.query_name
                if read_name is None:
                    continue

                nmatches = matches_from_cigar(read.cigartuples) if read.cigartuples else 0
                length = read.infer_read_length()
                if length is None:
                    length = len(read.seq) if read.seq else 0

                rec = RecordSummary(nmatches, length)
                file_reads.setdefault(read_name, {}) \
                    .setdefault(bam_file, []).append(rec)
                read_names_in_file.add(read_name)

            sam.close()

            self.files.append(bam_file)
            for read_name, file_records in file_reads.items():
                if read_name not in self.read_hash:
                    self.read_names.append(read_name)
                    self.read_names_set.add(read_name)
                self.read_hash[read_name] = file_records

            for rn in read_names_in_file:
                if rn not in self.read_names_set:
                    self.read_names.append(rn)
                    self.read_names_set.add(rn)

            self.read_count += 1

    def reset_reads(self):
        self.reads_to_keep = list(self.read_names)

    def filter_reads_to_keep(self):
        pop_size = len(self.files)
        reads_to_keep_temp = []
        for read_name in self.read_names:
            this_read_hash = self.read_hash.get(read_name, {})
            total_reads = len(this_read_hash)
            frequency = total_reads / pop_size if pop_size > 0 else 0.0
            if frequency >= self.frequency_threshold:
                reads_to_keep_temp.append(read_name)
        self.reads_to_keep = reads_to_keep_temp

    def sample_reads(self):
        rng = random.Random()

        if len(self.reads_to_keep) <= self.max_reads:
            return

        num_files = len(self.files)
        if self.max_reads <= num_files:
            self.reads_to_keep = rng.sample(self.reads_to_keep, self.max_reads)
            return

        min_per_file = self.max_reads // num_files
        selected = set()

        for file in self.files:
            candidates = [
                r for r in self.reads_to_keep
                if r not in selected
                and self.read_presence.get(r, {}).__contains__(file)
            ]
            take = min(min_per_file, len(candidates))
            chosen = rng.sample(candidates, take) if take > 0 else []
            for r in chosen:
                selected.add(r)

        remaining = self.max_reads - len(selected)
        if remaining > 0:
            pool = [r for r in self.reads_to_keep if r not in selected]
            chosen = rng.sample(pool, min(remaining, len(pool)))
            for r in chosen:
                selected.add(r)

        self.reads_to_keep = list(selected)

    def guard_file_coverage(self):
        rng = random.Random()
        for file in self.files:
            has_coverage = any(
                self.read_presence.get(r, {}).__contains__(file)
                for r in self.reads_to_keep
            )
            if has_coverage:
                continue

            candidates = [
                name for name, files in self.read_presence.items()
                if file in files and name not in self.reads_to_keep
            ]
            if candidates:
                read = rng.choice(candidates)
                self.reads_to_keep.append(read)
                self.read_names_set.add(read)
            else:
                print(f"[WARN] File '{file}' has zero reads in pool")

    def phase1_parse_file(self, bam_file, per_file_limit):
        file_count = 0
        try:
            sam = pysam.AlignmentFile(bam_file, "rb")
        except Exception as e:
            print(f"[WARN] Failed to open {bam_file}: {e}")
            return

        for read in sam:
            if read.flag == 4:
                continue
            read_name = read.query_name
            if read_name is None:
                continue

            if read_name not in self.read_presence:
                if file_count >= per_file_limit:
                    continue
                file_count += 1

            self.read_presence.setdefault(read_name, set()).add(bam_file)

        sam.close()

    def phase1_parse_files(self, bam_files, file_counts):
        total_mapped = sum(file_counts) if file_counts else 0

        for i, bam_file in enumerate(bam_files):
            self.files.append(bam_file)
            if total_mapped > 0:
                raw = (self.max_reads * file_counts[i]) // total_mapped
            else:
                raw = self.max_reads // len(bam_files)
            per_file_limit = max(1, raw)
            self.phase1_parse_file(bam_file, per_file_limit)
            self.phase1_file_counts[bam_file] = 0

    def phase2_parse_file(self, bam_file):
        reads_to_keep_set = set(self.reads_to_keep)
        try:
            sam = pysam.AlignmentFile(bam_file, "rb")
        except Exception as e:
            print(f"[WARN] Failed to open {bam_file}: {e}")
            return

        for read in sam:
            if read.flag == 4:
                continue
            read_name = read.query_name
            if read_name is None:
                continue
            if read_name not in reads_to_keep_set:
                continue

            nmatches = matches_from_cigar(read.cigartuples) if read.cigartuples else 0
            length = read.infer_read_length()
            if length is None:
                length = len(read.seq) if read.seq else 0

            rec = RecordSummary(nmatches, length)
            self.read_hash.setdefault(read_name, {}).setdefault(bam_file, []).append(rec)

        sam.close()

    def phase2_parse_files(self):
        for bam_file in self.files:
            self.phase2_parse_file(bam_file)

    def phase1_filter_reads_by_frequency(self):
        pop_size = len(self.files)
        self.reads_to_keep = [
            name for name, files in self.read_presence.items()
            if len(files) / pop_size >= self.frequency_threshold
        ]
        self.read_names = list(self.reads_to_keep)
        self.read_names_set = set(self.reads_to_keep)

    def parse_files(self, files, file_counts):
        self.phase1_parse_files(files, file_counts)

    def filter_samples_to_keep(self, keep_empty_samples=False):
        files_to_keep = []
        for sam_file in self.files:
            has_data = any(
                sam_file in v for v in self.read_hash.values()
            )
            if has_data:
                files_to_keep.append(sam_file)
            elif keep_empty_samples:
                print(f"[WARN] File '{sam_file}' has zero reads after Phase 2")
                files_to_keep.append(sam_file)
        self.files = files_to_keep

    def sample_total_reads_vector(self):
        file_count = len(self.files)
        counts = []
        for j in range(file_count):
            sam_file = self.files[j]
            cnt = sum(
                1 for read_name in self.reads_to_keep
                if sam_file in self.read_hash.get(read_name, {})
            )
            counts.append(cnt)
        return counts
