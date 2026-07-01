use ahash::{AHashMap, AHashSet};
use rand::seq::SliceRandom;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Reader, Record};
use std::cell::OnceCell;
use std::fs::File;
use std::io::prelude::*;

pub fn transpose_matrix(matrix: Vec<Vec<i32>>) -> Vec<Vec<i32>> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return vec![];
    }
    let num_cols = matrix[0].len();
    (0..num_cols)
        .into_par_iter()
        .map(|i| matrix.iter().map(|row| row[i]).collect())
        .collect()
}

pub fn transpose_matrix_float(matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return vec![];
    }
    let num_cols = matrix[0].len();
    (0..num_cols)
        .into_par_iter()
        .map(|i| matrix.iter().map(|row| row[i]).collect())
        .collect()
}

pub fn private_reads_rows(matrix1: &Vec<Vec<i32>>, matrix2: Vec<Vec<i32>>) -> i32 {
    let mut private_reads = 0;
    for i in 0..matrix1.len() {
        let sum1: i32 = matrix1[i].iter().sum();
        let sum2: i32 = matrix2[i].iter().sum();
        if sum1 - sum2 == 0 {
            private_reads += 1;
        }
    }
    private_reads
}

pub fn private_reads_proportion(
    subset_matrix1: &Vec<Vec<i32>>,
    matrix2: Vec<Vec<i32>>,
) -> (i32, f32) {
    let mut private_reads = 0;
    let mut matrix_1_reads = 0;

    for i in 0..subset_matrix1.len() {
        let subset_sum1: i32 = subset_matrix1[i].iter().sum();
        if subset_sum1 == 0 {
            continue;
        }
        let sum2: i32 = matrix2[i].iter().sum();
        if subset_sum1 - sum2 == 0 {
            private_reads += 1;
        }
        matrix_1_reads += 1;
    }

    let proportion = private_reads as f32 / matrix_1_reads as f32;
    (private_reads, proportion)
}

pub struct RecordEvalNumMatches {}

impl RecordEvalNumMatches {
    fn matches_from_cigar(cigar: &str) -> i32 {
        let mut matches = 0;
        let mut current_match = String::new();
        for c in cigar.chars() {
            if c.is_alphabetic() {
                if c == 'M' {
                    matches += current_match.parse::<i32>().unwrap();
                }
                current_match = String::new();
            } else {
                current_match.push(c);
            }
        }
        matches
    }

    fn eval_record(record: &Record) -> i32 {
        let cigar = record
            .cigar()
            .iter()
            .map(|x| x.to_string())
            .collect::<String>();
        RecordEvalNumMatches::matches_from_cigar(&cigar)
    }
}

pub struct RecordSummary {
    pub nmatches: i32,
    pub length: i32,
}

impl RecordSummary {
    pub fn percent_identity(&self) -> f32 {
        self.nmatches as f32 / self.length as f32
    }
}

pub struct MultipleSamParser {
    pub read_hash: AHashMap<String, AHashMap<String, Vec<RecordSummary>>>,
    pub _read_names: Vec<String>,
    pub _read_names_set: AHashSet<String>,
    pub files: Vec<String>,
    pub reads_to_keep: Vec<String>,
    pub read_count: usize,
    pub frequency_threshold: f32,
    pub max_reads: i32,
    pub _cached_total_reads: OnceCell<Vec<i32>>,
    pub _cached_presence_absence: OnceCell<Vec<Vec<i32>>>,
    pub read_presence: AHashMap<String, AHashSet<String>>,
    pub _phase1_file_counts: AHashMap<String, i32>,
}

pub struct CladeSummary {
    pub leaves: Vec<String>,
    pub name: String,
    pub private_mutations: i32,
    pub shared_mutations: i32,
    pub total_mutations: i32,
}

pub fn find_bam_files(directory: &str) -> Vec<String> {
    let mut bam_files: Vec<String> = Vec::new();
    let paths = std::fs::read_dir(directory).unwrap();

    for path in paths {
        let path = path.unwrap().path();
        let path_str = path.to_str().unwrap();
        if path_str.ends_with(".bam") {
            bam_files.push(path_str.to_string());
        }
    }

    bam_files
}

pub fn validate_bam_files(files: Vec<String>) -> Vec<String> {
    files
        .into_iter()
        .filter(|f| {
            let path = std::path::Path::new(f);
            if !path.exists() {
                println!("[WARN] File not found, skipping: {}", f);
                false
            } else if !f.ends_with(".bam") {
                println!("[WARN] Not a .bam file, skipping: {}", f);
                false
            } else {
                true
            }
        })
        .collect()
}

pub fn count_mapped_reads(files: &[String]) -> Vec<i32> {
    files
        .iter()
        .map(|file| {
            let mut bam = Reader::from_path(file).unwrap();
            let mut record = Record::new();
            let mut count = 0;
            while let Some(r) = bam.read(&mut record) {
                r.expect("Failed to parse record");
                if record.flags() != 4 {
                    count += 1;
                }
            }
            count
        })
        .collect()
}

fn parse_bam_file_rayon(
    bam_file: &str,
) -> (
    String,
    AHashMap<String, AHashMap<String, Vec<RecordSummary>>>,
    AHashSet<String>,
) {
    let bam_file_owned = bam_file.to_string();
    let mut bam = Reader::from_path(bam_file).unwrap();
    let mut record = Record::new();
    let mut read_data: AHashMap<String, AHashMap<String, Vec<RecordSummary>>> = AHashMap::new();
    let mut read_names: AHashSet<String> = AHashSet::new();

    while let Some(r) = bam.read(&mut record) {
        r.expect("Failed to parse record");

        let read_buf = record.qname();
        let read_name = String::from_utf8(read_buf.to_vec()).unwrap();

        let flag = record.flags();
        if flag == 4 {
            continue;
        }

        let record_summary = RecordSummary {
            nmatches: RecordEvalNumMatches::eval_record(&record),
            length: record.seq().len() as i32,
        };

        read_data
            .entry(read_name.clone())
            .or_insert_with(AHashMap::new)
            .entry(bam_file_owned.clone())
            .or_insert_with(Vec::new)
            .push(record_summary);

        read_names.insert(read_name);
    }

    (bam_file_owned, read_data, read_names)
}

impl MultipleSamParser {
    pub fn new(frequency_threshold: f32, max_reads: i32) -> MultipleSamParser {
        MultipleSamParser {
            read_hash: AHashMap::new(),
            _read_names: Vec::new(),
            _read_names_set: AHashSet::new(),
            files: Vec::new(),
            reads_to_keep: Vec::new(),
            read_count: 0,
            frequency_threshold,
            max_reads,
            _cached_total_reads: OnceCell::new(),
            _cached_presence_absence: OnceCell::new(),
            read_presence: AHashMap::new(),
            _phase1_file_counts: AHashMap::new(),
        }
    }

    pub fn parse_bam_files(&mut self, bam_files: Vec<String>) {
        let results: Vec<(
            String,
            AHashMap<String, AHashMap<String, Vec<RecordSummary>>>,
            AHashSet<String>,
        )> = bam_files
            .par_iter()
            .map(|f| parse_bam_file_rayon(f))
            .collect();

        for (file, read_data, read_names_set) in results {
            self.files.push(file);
            for (read_name, file_records) in read_data.into_iter() {
                if !self.read_hash.contains_key(&read_name) {
                    self._read_names.push(read_name.clone());
                    self._read_names_set.insert(read_name.clone());
                }
                self.read_hash.insert(read_name, file_records);
            }
            for rn in read_names_set {
                if !self._read_names_set.contains(&rn) {
                    self._read_names.push(rn.clone());
                    self._read_names_set.insert(rn);
                }
            }
            self.read_count += 1;
        }
    }

    pub fn parse_files(&mut self, files: Vec<String>, file_counts: Vec<i32>) {
        self.phase1_parse_files(files, file_counts);
    }

    fn phase1_parse_file(&mut self, bam_file: &str, per_file_limit: usize) {
        let mut bam = Reader::from_path(bam_file).unwrap();
        let mut record = Record::new();
        let bam_file_owned = bam_file.to_string();
        let mut file_count = 0;

        while let Some(r) = bam.read(&mut record) {
            r.expect("Failed to parse record");

            let flag = record.flags();
            if flag == 4 {
                continue; // skip unmapped reads
            }

            let read_buf = record.qname();
            let read_name = String::from_utf8(read_buf.to_vec()).unwrap();

            if !self.read_presence.contains_key(&read_name) {
                if file_count >= per_file_limit {
                    continue;
                }
                file_count += 1;
            }

            self.read_presence
                .entry(read_name)
                .or_insert_with(AHashSet::new)
                .insert(bam_file_owned.clone());
        }
    }

    pub fn phase1_parse_files(&mut self, bam_files: Vec<String>, file_counts: Vec<i32>) {
        let total_mapped: i64 = file_counts.iter().map(|&c| c as i64).sum();
        let max_reads = self.max_reads as i64;

        for (i, bam_file) in bam_files.iter().enumerate() {
            self.files.push(bam_file.clone());

            let raw = if total_mapped > 0 {
                (max_reads * file_counts[i] as i64) / total_mapped
            } else {
                max_reads / bam_files.len() as i64
            };
            let per_file_limit = std::cmp::max(1, raw as usize);

            self.phase1_parse_file(bam_file, per_file_limit);
            self._phase1_file_counts.insert(bam_file.clone(), 0);
        }
    }

    pub fn phase1_filter_reads_by_frequency(&mut self) {
        let pop_size = self.files.len() as f32;
        let threshold = self.frequency_threshold;

        self.reads_to_keep = self
            .read_presence
            .iter()
            .filter(|(_read_name, files)| {
                let frequency = files.len() as f32 / pop_size;
                frequency >= threshold
            })
            .map(|(name, _)| name.clone())
            .collect();

        self._read_names = self.reads_to_keep.clone();
        self._read_names_set = self.reads_to_keep.iter().cloned().collect();
    }

    pub fn phase2_parse_files(&mut self) {
        let files = self.files.clone();
        for bam_file in files {
            self.phase2_parse_file(&bam_file);
        }
    }

    fn phase2_parse_file(&mut self, bam_file: &str) {
        let mut bam = Reader::from_path(bam_file).unwrap();
        let mut record = Record::new();
        let bam_file_owned = bam_file.to_string();
        let reads_to_keep_set: AHashSet<&String> = self.reads_to_keep.iter().collect();

        while let Some(r) = bam.read(&mut record) {
            r.expect("Failed to parse record");

            let flag = record.flags();
            if flag == 4 {
                continue;
            }

            let read_buf = record.qname();
            let read_name = String::from_utf8(read_buf.to_vec()).unwrap();

            if !reads_to_keep_set.contains(&read_name) {
                continue;
            }

            let record_summary = RecordSummary {
                nmatches: RecordEvalNumMatches::eval_record(&record),
                length: record.seq().len() as i32,
            };

            self.read_hash
                .entry(read_name.clone())
                .or_insert_with(AHashMap::new)
                .entry(bam_file_owned.clone())
                .or_insert_with(Vec::new)
                .push(record_summary);
        }
    }

    pub fn reset_reads(&mut self) {
        self.reads_to_keep = self._read_names.clone();
    }

    pub fn sample_reads(&mut self) {
        let mut rng = rand::thread_rng();

        if self.reads_to_keep.len() <= self.max_reads as usize {
            return;
        }

        let num_files = self.files.len();
        if self.max_reads as usize <= num_files {
            self.reads_to_keep = self
                .reads_to_keep
                .choose_multiple(&mut rng, self.max_reads as usize)
                .cloned()
                .collect();
            return;
        }

        let min_per_file = self.max_reads as usize / num_files;
        let mut selected: AHashSet<String> = AHashSet::new();

        for file in &self.files {
            let candidates: Vec<&String> = self
                .reads_to_keep
                .iter()
                .filter(|r| !selected.contains(*r))
                .filter(|r| {
                    self.read_presence
                        .get(*r)
                        .map_or(false, |files| files.contains(file))
                })
                .collect();
            let take = std::cmp::min(min_per_file, candidates.len());
            for &r in candidates.choose_multiple(&mut rng, take) {
                selected.insert(r.clone());
            }
        }

        let remaining = self.max_reads as usize - selected.len();
        if remaining > 0 {
            let pool: Vec<&String> = self
                .reads_to_keep
                .iter()
                .filter(|r| !selected.contains(*r))
                .collect();
            for &r in pool.choose_multiple(&mut rng, remaining) {
                selected.insert(r.clone());
            }
        }

        self.reads_to_keep = selected.into_iter().collect();
    }

    pub fn guard_file_coverage(&mut self) {
        for file in &self.files {
            let has_coverage = self.reads_to_keep.iter().any(|r| {
                self.read_presence
                    .get(r)
                    .map_or(false, |files| files.contains(file))
            });
            if has_coverage {
                continue;
            }

            let candidates: Vec<&String> = self
                .read_presence
                .iter()
                .filter(|(_, files)| files.contains(file))
                .filter(|(name, _)| !self.reads_to_keep.contains(name))
                .map(|(name, _)| name)
                .collect();

            if let Some(read) = candidates.choose(&mut rand::thread_rng()) {
                println!(
                    "[INFO] Guard: injecting read '{}' to cover file '{}'",
                    read, file
                );
                self.reads_to_keep.push((*read).clone());
                self._read_names_set.insert((*read).clone());
            } else {
                println!(
                    "[WARN] File '{}' has zero reads in pool — cannot guarantee coverage",
                    file
                );
            }
        }
    }

    pub fn filter_reads_to_keep(&mut self) {
        let pop_size: f32 = self.files.len() as f32;
        let mut reads_to_keep_temp: Vec<String> = Vec::new();

        for read_name in &self._read_names {
            let this_read_hash = self.read_hash.get(read_name).unwrap();
            let total_reads = this_read_hash.len() as f32;
            let frequency = total_reads / pop_size;

            if frequency >= self.frequency_threshold {
                reads_to_keep_temp.push(read_name.clone());
            }
        }

        self.reads_to_keep = reads_to_keep_temp;
    }

    pub fn filter_samples_to_keep(&mut self, keep_empty_samples: bool) {
        let mut files_to_keep: Vec<String> = Vec::new();
        for sam_file in &self.files {
            if self.read_hash.values().any(|v| v.contains_key(sam_file)) {
                files_to_keep.push(sam_file.clone());
            } else {
                if keep_empty_samples {
                    println!(
                        "[WARN] File '{}' has zero reads after Phase 2 — including empty row in matrix",
                        sam_file
                    );
                    files_to_keep.push(sam_file.clone());
                }
            }
        }
        self.files = files_to_keep;
    }

    pub fn to_matrix_num_matches(&self) -> Vec<Vec<i32>> {
        let files = self.files.clone();
        let read_hash = &self.read_hash;
        let reads_to_keep = &self.reads_to_keep;

        reads_to_keep
            .par_iter()
            .map(|read_name| {
                let this_read_hash = read_hash.get(read_name).unwrap();
                let mut this_read_vec = Vec::with_capacity(files.len());
                for sam_file in &files {
                    if let Some(this_file_vec) = this_read_hash.get(sam_file) {
                        let max_value = this_file_vec.iter().map(|x| x.nmatches).max();
                        this_read_vec.push(max_value.unwrap_or(0));
                    } else {
                        this_read_vec.push(0);
                    }
                }
                this_read_vec
            })
            .collect()
    }

    pub fn to_matrix_presence_absence(&self) -> Vec<Vec<i32>> {
        if let Some(cached) = self._cached_presence_absence.get() {
            return cached.clone();
        }

        let files = self.files.clone();
        let read_hash = &self.read_hash;
        let reads_to_keep = &self.reads_to_keep;

        let result: Vec<Vec<i32>> = reads_to_keep
            .par_iter()
            .map(|read_name| {
                let this_read_hash = read_hash.get(read_name).unwrap();
                let mut this_read_vec = Vec::with_capacity(files.len());
                for sam_file in &files {
                    this_read_vec.push(if this_read_hash.contains_key(sam_file) {
                        1
                    } else {
                        0
                    });
                }
                this_read_vec
            })
            .collect();

        let _ = self._cached_presence_absence.set(result.clone());
        result
    }

    pub fn to_matrix_percent_match(&self) -> Vec<Vec<f32>> {
        let files = self.files.clone();
        let read_hash = &self.read_hash;
        let reads_to_keep = &self.reads_to_keep;

        reads_to_keep
            .par_iter()
            .map(|read_name| {
                let this_read_hash = read_hash.get(read_name).unwrap();
                let mut this_read_vec = Vec::with_capacity(files.len());
                for sam_file in &files {
                    if let Some(this_file_vec) = this_read_hash.get(sam_file) {
                        let max_value = this_file_vec
                            .iter()
                            .map(|x| x.percent_identity())
                            .fold(f32::NAN, |max, val| if val > max { val } else { max });
                        if max_value.is_nan() {
                            this_read_vec.push(0.0);
                        } else {
                            this_read_vec.push(max_value);
                        }
                    } else {
                        this_read_vec.push(0.0);
                    }
                }
                this_read_vec
            })
            .collect()
    }

    pub fn matrix_subset_files(&self, files: Vec<String>) -> Vec<Vec<i32>> {
        let read_hash = &self.read_hash;
        let reads_to_keep = &self.reads_to_keep;

        reads_to_keep
            .par_iter()
            .map(|read_name| {
                let this_read_hash = read_hash.get(read_name).unwrap();
                let mut this_read_vec = Vec::with_capacity(files.len());
                for sam_file in &files {
                    this_read_vec.push(if this_read_hash.contains_key(sam_file) {
                        1
                    } else {
                        0
                    });
                }
                this_read_vec
            })
            .collect()
    }

    pub fn subset_matrix_from_list(
        &self,
        full_reads_matrix: Vec<Vec<i32>>,
        names: Vec<String>,
        subset: Vec<String>,
    ) -> Vec<Vec<i32>> {
        let indices: Vec<usize> = subset
            .iter()
            .filter_map(|name| names.iter().position(|x| x == name))
            .collect();

        full_reads_matrix
            .par_iter()
            .map(|row| {
                let mut this_row = Vec::with_capacity(indices.len());
                for &index in &indices {
                    this_row.push(row[index]);
                }
                this_row
            })
            .collect()
    }

    pub fn standardize_matrix(&self, matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        matrix
            .par_iter()
            .map(|row| {
                let row_sum: f64 = row.iter().sum();
                if row_sum == 0.0 {
                    return vec![0.0; row.len()];
                }
                row.iter().map(|&val| val / row_sum).collect()
            })
            .collect()
    }

    pub fn sample_total_reads_vector(&self) -> Vec<i32> {
        if let Some(cached) = self._cached_total_reads.get() {
            return cached.clone();
        }

        let files = self.files.clone();
        let read_hash = &self.read_hash;
        let reads_to_keep = &self.reads_to_keep;
        let file_count = files.len();

        let counts: Vec<usize> = (0..file_count)
            .into_par_iter()
            .map(|j| {
                let sam_file = &files[j];
                reads_to_keep
                    .iter()
                    .filter(|read_name| {
                        if let Some(this_read_hash) = read_hash.get(*read_name) {
                            this_read_hash.contains_key(sam_file)
                        } else {
                            false
                        }
                    })
                    .count()
            })
            .collect();

        let total_reads_vector: Vec<i32> = counts.into_iter().map(|c| c as i32).collect();

        let _ = self._cached_total_reads.set(total_reads_vector.clone());
        total_reads_vector
    }

    pub fn standardize_by_files_total(&self, matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let total_reads_vector = self.sample_total_reads_vector();
        let num_cols = if matrix.is_empty() {
            0
        } else {
            matrix[0].len()
        };

        matrix
            .par_iter()
            .enumerate()
            .map(|(i, row)| {
                if total_reads_vector[i] == 0 {
                    return vec![0.0; num_cols];
                }
                row.par_iter()
                    .enumerate()
                    .map(|(j, &val)| {
                        if total_reads_vector[j] > 0 {
                            val / total_reads_vector[j] as f64
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect()
    }

    pub fn standardize_by_files_total_assymetric(&self, matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let total_reads_vector = self.sample_total_reads_vector();
        let num_cols = if matrix.is_empty() {
            0
        } else {
            matrix[0].len()
        };

        matrix
            .par_iter()
            .enumerate()
            .map(|(i, row)| {
                if total_reads_vector[i] == 0 {
                    return vec![0.0; num_cols];
                }
                row.par_iter()
                    .map(|&val| {
                        if total_reads_vector[i] > 0 {
                            val / total_reads_vector[i] as f64
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect()
    }

    pub fn invert_matrix(&self, matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        matrix
            .par_iter()
            .map(|row| row.iter().map(|&val| 1.0 - val).collect())
            .collect()
    }

    pub fn to_distance_matrix_int(&self, matrix: Vec<Vec<i32>>) -> Vec<Vec<f64>> {
        let transposed_matrix = transpose_matrix(matrix);
        let n = transposed_matrix.len();
        let m = if n > 0 { transposed_matrix[0].len() } else { 0 };

        (0..n)
            .into_par_iter()
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let mut distance = 0;
                        for k in 0..m {
                            distance += (transposed_matrix[i][k] - transposed_matrix[j][k]).abs();
                        }
                        distance as f64
                    })
                    .collect()
            })
            .collect()
    }

    pub fn to_distance_matrix_float(&self, matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let transposed_matrix = transpose_matrix_float(matrix);
        let n = transposed_matrix.len();
        let m = if n > 0 { transposed_matrix[0].len() } else { 0 };

        (0..n)
            .into_par_iter()
            .map(|i| {
                (0..n)
                    .map(|j| {
                        let mut distance = 0.0;
                        for k in 0..m {
                            distance += (transposed_matrix[i][k] - transposed_matrix[j][k]).abs();
                        }
                        distance
                    })
                    .collect()
            })
            .collect()
    }

    pub fn to_shared_mut_matrix(&self) -> Vec<Vec<f64>> {
        let files = &self.files;
        let file_count = files.len();
        let mut matrix = vec![vec![0.0_f64; file_count]; file_count];

        for read_name in &self.reads_to_keep {
            if let Some(this_read_hash) = self.read_hash.get(read_name) {
                let present_files: Vec<usize> = files
                    .iter()
                    .enumerate()
                    .filter(|(_, f)| this_read_hash.contains_key(*f))
                    .map(|(i, _)| i)
                    .collect();

                for i in 0..present_files.len() {
                    for j in i..present_files.len() {
                        matrix[present_files[i]][present_files[j]] += 1.0;
                        if i != j {
                            matrix[present_files[j]][present_files[i]] += 1.0;
                        }
                    }
                }
            }
        }

        matrix
    }

    pub fn to_file_reads(
        &self,
        file_name: &str,
        matrix: &Vec<Vec<i32>>,
        sep: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut file = File::create(file_name)?;

        let header = self.files.join(sep);
        file.write(String::from("read_name").as_bytes())?;
        file.write(sep.as_bytes())?;
        file.write(header.as_bytes())?;
        file.write_all(b"\n")?;

        for (i, row) in matrix.iter().enumerate() {
            let read_name = self.reads_to_keep[i].clone();
            let row_str = row
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(sep);
            file.write(read_name.as_bytes())?;
            file.write_all(sep.as_bytes())?;
            file.write_all(row_str.as_bytes())?;
            file.write_all(b"\n")?;
        }

        Ok(())
    }
}
