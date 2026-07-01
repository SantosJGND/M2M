pub mod core;
pub mod phylo;
pub mod read_phylo;
pub mod sam_utils;
use core::cli::Cli;
use phylo::phylogeny::neighbor_joining;
use read_phylo::read_phylogeny::{
    directory_exists_create, output_reports, print_distance_matrix, print_presence_absence_matrix,
    print_tree_edges, print_tree_newick, recursive_get_private_nodes,
    write_all_node_statistics,
};
use sam_utils::sam_parser::{
    count_mapped_reads, find_bam_files, transpose_matrix, validate_bam_files, MultipleSamParser,
};
use speedytree::DistanceMatrix;

fn main() {
    let config = Cli::new();
    let mut bam_files = config._files_list;
    let max_reads = config.max_reads;
    let _keep_empty_samples = config._keep_empty_samples;
    let input_directory = config.input_directory;
    let frequency_threshold = config.frequency_threshold;
    let _adjusted_distances = config.adjusted_distances;
    let private_threshold = config.threshold;
    let output_directory = config.output_directory;
    let cluster_analysis = config.cluster_analysis;

    if bam_files.len() > 0 && input_directory.len() > 0 {
        println!("Please provide either a list of BAM files or an input directory, not both.");
        return;
    }
    if bam_files.len() == 0 && input_directory.len() == 0 {
        println!("Please provide a list of BAM files or an input directory.");
        return;
    }

    if bam_files.len() > 0 {
        println!("Using BAM files provided in the command line arguments.");
    } else {
        bam_files = find_bam_files(&input_directory);

        if bam_files.len() == 0 {
            println!("No bam files found in directory. Exiting.");
            return;
        }
    }

    bam_files = validate_bam_files(bam_files);
    if bam_files.len() == 0 {
        println!("No valid BAM files found. Exiting.");
        return;
    }

    let mut parser = MultipleSamParser::new(frequency_threshold, max_reads);

    let file_counts = count_mapped_reads(&bam_files);
    println!("[INFO] Mapped read counts per file: {:?}", file_counts);

    parser.parse_files(bam_files, file_counts);
    println!("[INFO] BAM files parsed (Phase 1).");
    parser.phase1_filter_reads_by_frequency();
    println!("[INFO] Reads filtered by frequency.");
    parser.sample_reads();
    println!("[INFO] Reads sampled (limiting to max_reads for Phase 2).");
    parser.guard_file_coverage();
    println!("[INFO] File coverage guarded.");
    parser.phase2_parse_files();
    println!("[INFO] BAM files parsed (Phase 2 - sampled reads only).");
    parser.filter_samples_to_keep(_keep_empty_samples);
    println!("[INFO] Samples filtered.");
    if parser.reads_to_keep.len() == 0 {
        println!("No reads to process after filtering. Exiting.");
        return;
    }

    let matrix = parser.to_shared_mut_matrix();
    println!(
        "[INFO] Distance matrix built with {} samples.",
        matrix.len()
    );

    // break with warning if distance matrix is empty
    if matrix.len() == 0 {
        println!("Distance matrix is empty. Exiting.");
        return;
    }

    if matrix[0].len() == 0 {
        println!("Distance matrix is empty. Exiting.");
        return;
    }

    let shared_matrix: Vec<Vec<f64>> = parser.standardize_by_files_total(&matrix);
    let inverted_distance_matrix = parser.invert_matrix(shared_matrix.clone());

    directory_exists_create(&output_directory).unwrap();

    print_distance_matrix(
        &inverted_distance_matrix,
        &parser.files,
        output_directory.clone(),
    );

    let presence_absence_matrix = parser.to_matrix_presence_absence();
    println!(
        "[INFO] Presence/absence matrix built with {} samples and {} nodes.",
        presence_absence_matrix.len(),
        presence_absence_matrix[0].len()
    );

    let transposed_pa_matrix = transpose_matrix(presence_absence_matrix.clone());
    print_presence_absence_matrix(
        &transposed_pa_matrix,
        &parser.files,
        output_directory.clone(),
    );

    if !cluster_analysis {
        println!("[INFO] Cluster analysis disabled. Skipping tree and clade reports.");
        return;
    }

    let shared_matrix_assymetric = parser.standardize_by_files_total_assymetric(matrix);

    let phylo = neighbor_joining(
        DistanceMatrix::build(inverted_distance_matrix.clone(), parser.files.clone()).unwrap(),
    );

    println!(
        "[INFO] NJ tree built with {} leaves and {} nodes.",
        shared_matrix.len(),
        phylo.get_all_nodes().len()
    );

    let private_nodes = recursive_get_private_nodes(
        phylo.clone(),
        &parser,
        &private_threshold,
        presence_absence_matrix.clone(),
        DistanceMatrix {
            matrix: shared_matrix_assymetric.clone(),
            names: parser.files.clone(),
        },
        parser.files.clone(),
    );

    print_tree_edges(&phylo, output_directory.clone());
    print_tree_newick(&phylo, output_directory.clone());

    _ = write_all_node_statistics(
        phylo,
        &parser,
        presence_absence_matrix,
        DistanceMatrix {
            matrix: shared_matrix_assymetric.clone(),
            names: parser.files.clone(),
        },
        parser.files.clone(),
        output_directory.clone(),
    );

    _ = output_reports(private_nodes, output_directory, parser.files.clone());
    println!("[INFO] Reports written to output directory.");
}
