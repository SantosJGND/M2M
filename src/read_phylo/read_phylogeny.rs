use crate::phylo::phylogeny::{min_pairwise_shared, Phylogeny};
use crate::sam_utils::sam_parser::{private_reads_proportion, MultipleSamParser};
use speedytree::DistanceMatrix;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
#[derive(Debug, Clone)]
pub struct PhyloPrivateProportion {
    pub node: Phylogeny,
    pub private_proportion: f32,
    pub private_reads: i32,
    pub min_pairwise_dist: f64,
}

// recursive function to get nodes with private proportion > 0.1.
pub fn get_phylo_private_proportion(
    node: Phylogeny,
    parser: &MultipleSamParser,
    full_reads_matrix: Vec<Vec<i32>>,
    names: Vec<String>,
) -> (i32, f32) {
    //
    //let full_reads_matrix = parser.to_matrix_presence_absence();
    let leaves = node.recursive_get_leaves_string();
    let leaves_list: Vec<String> = leaves.iter().cloned().collect();

    //let subset_matrix = parser.matrix_subset_files(leaves_list);
    let subset_matrix =
        parser.subset_matrix_from_list(full_reads_matrix.clone(), names, leaves_list);

    let private_data = private_reads_proportion(&subset_matrix, full_reads_matrix.clone());
    private_data
}

// all node statistics
pub fn all_node_statistics(
    node: Phylogeny,
    parser: &MultipleSamParser,
    full_reads_matrix: Vec<Vec<i32>>,
    distance_matrix: DistanceMatrix,
    names: Vec<String>,
) -> std::io::Result<Vec<(String, i32, f32, f64, usize)>> {
    let mut node_stats_vec = Vec::new();

    let leaves = node.recursive_get_leaves_string();
    let leaves_list: Vec<String> = leaves.iter().cloned().collect();

    let private_data = get_phylo_private_proportion(
        node.clone(),
        parser,
        full_reads_matrix.clone(),
        names.clone(),
    );
    let min_pairwise_dist = min_pairwise_shared(
        &distance_matrix,
        &node.recursive_get_leaves_string().into_iter().collect(),
    );
    node_stats_vec.push((
        node.name.clone(),
        private_data.0,
        private_data.1,
        min_pairwise_dist,
        leaves_list.len(),
    ));
    for child in node.children.iter() {
        let mut child_stats = all_node_statistics(
            child.0.clone(),
            parser,
            full_reads_matrix.clone(),
            distance_matrix.clone(),
            names.clone(),
        )?;
        node_stats_vec.append(&mut child_stats);
    }

    Ok(node_stats_vec)
}

pub fn write_all_node_statistics(
    node: Phylogeny,
    parser: &MultipleSamParser,
    full_reads_matrix: Vec<Vec<i32>>,
    distance_matrix: DistanceMatrix,
    names: Vec<String>,
    output_dir: String,
) -> std::io::Result<()> {
    let all_node_stats =
        all_node_statistics(node, parser, full_reads_matrix, distance_matrix, names)?;

    let mut file = File::create(format!("{}/all_node_statistics.tsv", output_dir).to_string())?;
    writeln!(
        file,
        "Node\tPrivate_Reads\tPrivate_Proportion\tMin_Pairwise_Dist\tNum_Leaves"
    )?;
    for (node_name, private_reads, private_proportion, min_pairwise_dist, num_leaves) in
        all_node_stats
    {
        writeln!(
            file,
            "{}\t{}\t{:.3}\t{:.3}\t{}",
            node_name, private_reads, private_proportion, min_pairwise_dist, num_leaves
        )?;
    }

    Ok(())
}

// recursive function to get nodes with private proportion > threshold.
pub fn recursive_get_private_nodes(
    node: Phylogeny,
    parser: &MultipleSamParser,
    threshold: &f64,
    full_reads_matrix: Vec<Vec<i32>>,
    distance_matrix: DistanceMatrix,
    names: Vec<String>,
) -> Vec<PhyloPrivateProportion> {
    let mut private_nodes: Vec<PhyloPrivateProportion> = Vec::new();

    if node.children.len() == 0 {
        let private_data = get_phylo_private_proportion(
            node.clone(),
            parser,
            full_reads_matrix.clone(),
            names.clone(),
        );

        let min_pairwise_dist = min_pairwise_shared(
            &distance_matrix.clone(),
            &node.recursive_get_leaves_string().into_iter().collect(),
        );

        // check if node not already in list:
        if !private_nodes.iter().any(|n| n.node.name == node.name) {
            private_nodes.push(PhyloPrivateProportion {
                node: node.clone(),
                private_proportion: private_data.1,
                private_reads: private_data.0,
                min_pairwise_dist,
            });
        }
    } else {
        for child in node.children.iter() {
            let private_data = get_phylo_private_proportion(
                child.0.clone(),
                parser,
                full_reads_matrix.clone(),
                names.clone(),
            );

            let min_pairwise_dist = min_pairwise_shared(
                &distance_matrix.clone(),
                &child.0.recursive_get_leaves_string().into_iter().collect(),
            );

            if min_pairwise_dist >= *threshold {
                // check if node not already in list:
                private_nodes.push(PhyloPrivateProportion {
                    node: child.0.clone(),
                    private_proportion: private_data.1,
                    private_reads: private_data.0,
                    min_pairwise_dist,
                });
                println!(
                    "[INFO] Node: {} added with private_proportion = {:.3}",
                    child.0.name, private_data.1
                );
            } else {
                let child_private_nodes = recursive_get_private_nodes(
                    child.0.clone(),
                    parser,
                    threshold,
                    full_reads_matrix.clone(),
                    distance_matrix.clone(),
                    names.clone(),
                );
                for child_node in child_private_nodes {
                    private_nodes.push(child_node);
                }
            }
        }
    }
    private_nodes
}

// Build a map from node name to parent node reference
fn build_parent_map<'a>(node: &'a Phylogeny, map: &mut HashMap<String, &'a Phylogeny>) {
    for (child, _) in &node.children {
        map.insert(child.name.clone(), node);
        build_parent_map(child, map);
    }
}

pub fn print_tree_edges(phylo: &Phylogeny, output_dir: String) {
    let mut parent_child_pairs = HashMap::new();
    build_parent_map(phylo, &mut parent_child_pairs);

    let mut file = File::create(format!("{}/nj_tree_edges.txt", output_dir).to_string()).unwrap();
    for (child, parent) in parent_child_pairs {
        writeln!(file, "{}\t{}", parent.name, child).unwrap();
    }
}

pub fn traverse_from_leaves_to_root(
    tree: &Phylogeny,
    parser: &MultipleSamParser,
    threshold: f64,
    full_reads_matrix: Vec<Vec<i32>>,
    distance_matrix: DistanceMatrix,
    names: Vec<String>,
) -> Vec<PhyloPrivateProportion> {
    let mut private_nodes: Vec<PhyloPrivateProportion> = Vec::new();

    // Build parent map
    let mut parent_map: HashMap<String, &Phylogeny> = HashMap::new();
    build_parent_map(tree, &mut parent_map);

    // Get all leaves of the tree
    let leaves = tree.recursive_get_leaves();
    for leaf in leaves {
        let mut current_node = &leaf;

        // Traverse upward until the root or threshold is exceeded

        loop {
            let private_data = get_phylo_private_proportion(
                current_node.clone(),
                parser,
                full_reads_matrix.clone(),
                names.clone(),
            );

            let min_pairwise_dist = min_pairwise_shared(
                &distance_matrix.clone(),
                &current_node
                    .recursive_get_leaves_string()
                    .into_iter()
                    .collect(),
            );

            //let min_pairwise_dist = min_pairwise_distance(
            //    &distance_matrix.clone(),
            //    &current_node
            //        .recursive_get_leaves_string()
            //        .into_iter()
            //        .collect(),
            //);

            if private_data.1 > threshold as f32 {
                // check if all children in encompassed children
                if !private_nodes
                    .iter()
                    .any(|n| n.node.name == current_node.name)
                {
                    private_nodes.push(PhyloPrivateProportion {
                        node: current_node.clone(),
                        private_proportion: private_data.1,
                        private_reads: private_data.0,
                        min_pairwise_dist,
                    });
                }

                break;

                //break;
            }

            // Move to parent using the map
            match parent_map.get(&current_node.name) {
                Some(parent) => {
                    current_node = parent;
                }
                None => {
                    break;
                }
            }
        }
    }
    let mut private_nodes_clone: Vec<PhyloPrivateProportion> = Vec::new();
    for node in private_nodes.iter() {
        let mut remove = false;
        let leaves = node.node.recursive_get_leaves_string();

        for onode in private_nodes.iter() {
            if node.node.name != onode.node.name {
                let olveaves = onode.node.recursive_get_leaves_string();

                if leaves.is_superset(&olveaves) {
                    remove = true;
                }
                use std::collections::HashSet;

                if leaves.len() != olveaves.len() {
                    continue;
                }

                // Assuming leaves1 and leaves2 are Vec<String>
                let set1: HashSet<_> = leaves.iter().map(|name| name.to_string()).collect();
                let set2: HashSet<_> = olveaves.iter().map(|name| name.to_string()).collect();

                if set1 == set2 {
                    remove = true;
                }
            }
        }

        if !remove {
            private_nodes_clone.push(node.clone());
        }
    }

    private_nodes_clone
}

pub fn directory_exists_create(directory: &str) -> std::io::Result<()> {
    if !std::path::Path::new(directory).exists() {
        std::fs::create_dir(directory)?;
    }
    Ok(())
}

// print distance matrix
pub fn print_distance_matrix(matrix: &Vec<Vec<f64>>, names: &Vec<String>, output_dir: String) {
    let mut file = File::create(format!("{}/distance_matrix.tsv", output_dir).to_string()).unwrap();
    writeln!(file, "\t{}", names.join("\t")).unwrap();
    for (i, row) in matrix.iter().enumerate() {
        writeln!(
            file,
            "{}\t{}",
            names[i],
            row.iter()
                .map(|x| format!("{:.3}", x))
                .collect::<Vec<_>>()
                .join("\t")
        )
        .unwrap();
    }
}

pub fn print_tree_newick(phylo: &Phylogeny, output_dir: String) {
    let mut file = File::create(format!("{}/nj_tree.newick", output_dir).to_string()).unwrap();
    writeln!(file, "{}", phylo.to_newick()).unwrap();
}

/*
* Print the presence/absence matrix
*/
pub fn print_presence_absence_matrix(
    matrix: &Vec<Vec<i32>>,
    names: &Vec<String>,
    output_dir: String,
) {
    let mut file =
        File::create(format!("{}/presence_absence_matrix.tsv", output_dir).to_string()).unwrap();
    writeln!(file, "\t{}", names.join("\t")).unwrap();
    for (i, row) in matrix.iter().enumerate() {
        writeln!(
            file,
            "{}\t{}",
            names[i],
            row.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join("\t")
        )
        .unwrap();
    }
}

// output reports
// one report for samples, two columns: sample name, node_name
// one report for clades, three columns: clade name, private reads, private proportion, comma separated list of leaves
pub fn output_reports(
    private_nodes: Vec<PhyloPrivateProportion>,
    output_dir: String,
    samples: Vec<String>,
) -> std::io::Result<()> {
    let mut sample_report = File::create(format!("{}/sample_report.tsv", output_dir).to_string())?;
    let mut clade_report = File::create(format!("{}/clade_report.tsv", output_dir).to_string())?;
    let node_name_prefix = "clade_";
    let mut private_nodes = private_nodes;
    if private_nodes.len() == 0 {
        return Ok(());
    }
    let mut samples_clustered: Vec<String> = Vec::new();

    for (i, node) in private_nodes.iter_mut().enumerate() {
        let new_name = format!("{}{}", node_name_prefix, i);
        let leaves = node.node.recursive_get_leaves_string();

        let leaves_list: Vec<String> = leaves.iter().cloned().collect();

        for leaf in &leaves_list {
            let sample_row_string = format!("{}\t{}", leaf, new_name);
            writeln!(sample_report, "{}", sample_row_string)?;

            samples_clustered.push(leaf.to_string());
        }

        writeln!(
            clade_report,
            "{}\t{}\t{}\t{}\t{}\t{}",
            new_name,
            node.private_reads,
            node.private_proportion,
            node.min_pairwise_dist,
            leaves_list.len(),
            leaves_list.join(",")
        )?;
    }

    for sample in samples {
        if !samples_clustered.contains(&sample) {
            writeln!(sample_report, "{}\t{}", sample, "unclustered")?;
        }
    }

    Ok(())
}
