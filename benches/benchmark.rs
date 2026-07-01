use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mapping_to_matrix::sam_utils::sam_parser::{find_bam_files, MultipleSamParser};

fn criterion_benchmark(c: &mut Criterion) {
    let test_data_dir = "test_data_larger";

    c.bench_function("parse_bam_files", |b| {
        b.iter(|| {
            let bam_files = find_bam_files(black_box(test_data_dir));
            let mut parser = MultipleSamParser::new(0.1, 500000);
            parser.parse_bam_files(black_box(bam_files));
            parser.reset_reads();
            parser.filter_reads_to_keep();
            parser.sample_reads();
            parser.filter_samples_to_keep(false);
            let read_count = parser.read_count;
            let files_count = parser.files.len();
            black_box((read_count, files_count))
        });
    });

    c.bench_function("to_matrix_presence_absence", |b| {
        b.iter(|| {
            let bam_files = find_bam_files(test_data_dir);
            let mut parser = MultipleSamParser::new(0.1, 500000);
            parser.parse_bam_files(bam_files);
            parser.reset_reads();
            parser.filter_reads_to_keep();
            parser.sample_reads();
            parser.filter_samples_to_keep(false);
            let matrix = parser.to_matrix_presence_absence();
            let dims = (matrix.len(), matrix.first().map(|r| r.len()).unwrap_or(0));
            black_box(dims)
        });
    });

    c.bench_function("to_shared_mut_matrix", |b| {
        b.iter(|| {
            let bam_files = find_bam_files(test_data_dir);
            let mut parser = MultipleSamParser::new(0.1, 500000);
            parser.parse_bam_files(bam_files);
            parser.reset_reads();
            parser.filter_reads_to_keep();
            parser.sample_reads();
            parser.filter_samples_to_keep(false);
            let matrix = parser.to_shared_mut_matrix();
            let dims = (matrix.len(), matrix.first().map(|r| r.len()).unwrap_or(0));
            black_box(dims)
        });
    });

    c.bench_function("standardize_by_files_total", |b| {
        b.iter(|| {
            let bam_files = find_bam_files(test_data_dir);
            let mut parser = MultipleSamParser::new(0.1, 500000);
            parser.parse_bam_files(bam_files);
            parser.reset_reads();
            parser.filter_reads_to_keep();
            parser.sample_reads();
            parser.filter_samples_to_keep(false);
            let matrix = parser.to_shared_mut_matrix();
            let standardized = parser.standardize_by_files_total(&matrix);
            let dims = (
                standardized.len(),
                standardized.first().map(|r| r.len()).unwrap_or(0),
            );
            black_box(dims)
        });
    });

    c.bench_function("full_pipeline", |b| {
        b.iter(|| {
            let bam_files = find_bam_files(test_data_dir);
            let mut parser = MultipleSamParser::new(0.1, 500000);
            parser.parse_bam_files(bam_files);
            parser.reset_reads();
            parser.filter_reads_to_keep();
            parser.sample_reads();
            parser.filter_samples_to_keep(false);

            if parser.reads_to_keep.is_empty() {
                return;
            }

            let matrix = parser.to_shared_mut_matrix();
            let shared_matrix = parser.standardize_by_files_total(&matrix);
            let inverted_distance_matrix = parser.invert_matrix(shared_matrix);

            let dims = (
                inverted_distance_matrix.len(),
                inverted_distance_matrix
                    .first()
                    .map(|r| r.len())
                    .unwrap_or(0),
            );
            black_box(dims);
        });
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
