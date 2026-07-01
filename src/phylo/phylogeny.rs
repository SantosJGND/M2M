use ordered_float::NotNan;
use speedytree::DistanceMatrix;
use std::cmp::{max, min};
use std::collections::HashSet;

use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct Phylogeny {
    /// The name of the current node.
    ///
    /// Can be empty for internal nodes.
    pub name: String,

    /// The children of the current node.
    ///
    /// Empty for leafs, and distances to the parent are optional.
    pub children: Vec<(Phylogeny, Option<NotNan<f64>>)>,
}

impl Phylogeny {
    pub fn to_newick(&self) -> String {
        if self.children.is_empty() {
            return self.name.clone();
        }

        let children_str: Vec<String> = self
            .children
            .iter()
            .map(|(child, dist)| {
                if let Some(d) = dist {
                    format!("{}:{}", child.to_newick(), d)
                } else {
                    child.to_newick()
                }
            })
            .collect();

        format!("({}){}", children_str.join(","), self.name)
    }

    /// Create a new phylogeny node from a name and list of children.
    pub fn new(name: &str, children: Vec<(Phylogeny, Option<NotNan<f64>>)>) -> Phylogeny {
        let phylo = Phylogeny {
            name: name.to_string(),
            children,
        };

        phylo
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }

    /// Create a new leaf node.
    pub fn new_leaf(name: &str) -> Phylogeny {
        Phylogeny::new(name, vec![])
    }

    pub fn clone(&self) -> Phylogeny {
        Phylogeny {
            name: self.name.clone(),
            children: self
                .children
                .iter()
                .map(|(child, dist)| (child.clone(), *dist))
                .collect(),
        }
    }

    /// Join two trees into a new named parent node, with given edge lengths.
    pub fn join_with_name(name: &str, l: Phylogeny, dl: f64, r: Phylogeny, dr: f64) -> Phylogeny {
        Phylogeny::new(
            name,
            vec![
                (l, Some(NotNan::new(dl).unwrap())),
                (r, Some(NotNan::new(dr).unwrap())),
            ],
        )
    }

    /// Join two trees into a new anonymous parent node, with given edge lengths.
    pub fn join(l: Phylogeny, dl: f64, r: Phylogeny, dr: f64, name: &str) -> Phylogeny {
        Phylogeny::join_with_name(name, l, dl, r, dr)
    }

    pub fn recursive_get_leaves_string(&self) -> HashSet<String> {
        let mut leaves: HashSet<String> = HashSet::new();

        if self.children.len() == 0 {
            leaves.insert(self.name.clone());
        } else {
            for child in self.children.iter() {
                let child_leaves = child.0.recursive_get_leaves_string();
                leaves.extend(child_leaves);
            }
        }

        leaves
    }
    pub fn recursive_get_leaves(&self) -> HashSet<Phylogeny> {
        let mut leaves: HashSet<Phylogeny> = HashSet::new();

        if self.children.len() == 0 {
            leaves.insert(self.clone());
        } else {
            for child in self.children.iter() {
                let child_leaves = child.0.recursive_get_leaves();
                leaves.extend(child_leaves.into_iter()); // Fix: Use `.into_iter()` here
            }
        }

        leaves
    }
    pub fn get_all_nodes(&self) -> Vec<Phylogeny> {
        let mut nodes: Vec<Phylogeny> = Vec::new();
        nodes.push(self.clone());

        for child in self.children.iter() {
            nodes.extend(child.0.get_all_nodes());
        }

        nodes
    }
}

pub fn min_pairwise_shared(distance_matrix: &DistanceMatrix, subset_names: &Vec<String>) -> f64 {
    let matrix_names = &distance_matrix.names;
    let matrix = &distance_matrix.matrix;

    let indices: Vec<usize> = subset_names
        .iter()
        .filter_map(|name| matrix_names.iter().position(|n| n == name))
        .collect();

    if indices.is_empty() {
        return f64::INFINITY;
    }

    if indices.len() == 1 {
        return 0.0;
    }

    let mut min_distance = f64::INFINITY;

    for i in 0..indices.len() {
        for j in i + 1..indices.len() {
            let dist = matrix[indices[i]][indices[j]]
                .min(matrix[indices[j]][indices[i]]);
            if dist < min_distance {
                min_distance = dist;
            }
        }
    }

    min_distance
}

pub fn neighbor_joining(distances_matrix: DistanceMatrix) -> Phylogeny {
    // Tuples of (phylogeny, number of leafs in tree, distance to leaf delta)
    let mut parts: Vec<Option<Phylogeny>> = distances_matrix
        .names
        .iter()
        .map(|name| Some(Phylogeny::new_leaf(name)))
        .collect();

    let mut distances = distances_matrix.matrix.clone();

    let mut inner_nodes_count = 0;

    while parts.iter().filter(|x| x.is_some()).count() > 2 {
        // Collect indices of non-None parts
        let indices: Vec<usize> = parts
            .iter()
            .enumerate()
            .filter_map(|(idx, x)| if x.is_some() { Some(idx) } else { None })
            .collect();

        let n = indices.len();

        //let sum_d = |i: usize| -> f64 {
        //    indices.iter().map(|&k| distances[i][k]).sum::<f64>() / (n as f64 - 2.)
        //};
        let sum_d = |i: usize| -> f64 {
            if n > 2 {
                indices.iter().map(|&k| distances[i][k]).sum::<f64>() / (n as f64 - 2.)
            } else {
                0.0 // Default to 0.0 if n <= 2
            }
        };

        use std::cmp::min;

        // Find pair (i, j) to merge
        let mut min_q = None;
        let mut min_pair = (0, 0);

        for &i in &indices {
            for &j in &indices {
                if i == j {
                    continue;
                }

                let mut min_here = Some(distances[i][j]);
                if distances[j][i] < distances[i][j] {
                    min_here = Some(distances[j][i]);
                }
                if min_q.is_none() || (min_here.is_some() && min_here < min_q) {
                    min_q = min_here;
                    min_pair = (i, j);
                }
            }
        }
        let (i, j) = (min(min_pair.0, min_pair.1), max(min_pair.0, min_pair.1));

        let di = distances[i][j] / 2. + (sum_d(i) - sum_d(j)) / (2. * (n as f64 - 2.));
        let di = if di.is_nan() { 0.0 } else { di };

        let dj = distances[i][j] - di;
        let dj = if dj.is_nan() { 0.0 } else { dj };

        let inner_node_name = format!("inner_{}", inner_nodes_count);
        let parent = Phylogeny::join(
            parts[i].take().unwrap(),
            di,
            parts[j].take().unwrap(),
            dj,
            &inner_node_name,
        );
        parts.push(Some(parent));
        inner_nodes_count += 1;

        // Build new row for the new node
        let mut new_row: Vec<f64> = Vec::with_capacity(distances.len() + 1);
        for &k in &indices {
            if k != i && k != j {
                let dk = (distances[i][k] + distances[j][k] - distances[i][j]) / 2.;
                new_row.push(dk);
            }
        }
        new_row.push(0.0);
        let mut row_use: Vec<f64> = new_row.clone();

        // Remove i and j by setting them to None
        parts[i] = None;
        parts[j] = None;

        // Update distances matrix
        let mut new_distances = Vec::new();
        for (idx, row) in distances.iter().enumerate() {
            if parts[idx].is_some() {
                let mut updated_row = Vec::new();
                for (jdx, &val) in row.iter().enumerate() {
                    if parts[jdx].is_some() {
                        updated_row.push(val);
                    }
                }
                new_distances.push(updated_row);
            }
        }
        for row in new_distances.iter_mut() {
            row.push(row_use.remove(0));
        }
        new_distances.push(new_row);

        distances = new_distances;
        // remove None from parts
        parts.retain(|x| x.is_some());
    }
    // Merge the two remaining vertices
    let indices: Vec<usize> = parts
        .iter()
        .enumerate()
        .filter_map(|(idx, x)| if x.is_some() { Some(idx) } else { None })
        .collect();
    if indices.len() == 2 {
        let inner_node_name = format!("inner_{}", inner_nodes_count);
        let i = indices[0];
        let j = indices[1];
        let d = distances[i][j] / 2.;
        parts[i] = Some(Phylogeny::join(
            parts[i].take().unwrap(),
            d,
            parts[j].take().unwrap(),
            d,
            &inner_node_name,
        ));
        parts[j] = None;
    }

    // Return the only remaining non-None tree
    parts
        .into_iter()
        .find_map(|x| x)
        .expect("Final tree should not be None")
}

pub fn upgma(distances: DistanceMatrix) -> Phylogeny {
    // Tuples of (phylogeny, number of leafs in tree, distance to leaf)
    let mut parts: Vec<Option<(Phylogeny, usize, f64)>> = distances
        .names
        .iter()
        .map(|name| Some((Phylogeny::new_leaf(name), 1, 0.0)))
        .collect();
    let mut distances: Vec<Vec<f64>> = distances.matrix.clone();
    let n = distances.len();

    let mut heap = BinaryHeap::with_capacity(n * n);
    for (i, ds) in distances.iter().enumerate() {
        for (j, &d) in ds.iter().enumerate() {
            if i == j {
                continue;
            }
            heap.push(Reverse((NotNan::new(d).unwrap(), i, j)));
        }
    }

    let mut inner_nodes_count = 0;

    while let Some(Reverse((d, i, j))) = heap.pop() {
        assert!(i != j);
        if parts[i].is_none() || parts[j].is_none() || distances[i][j] != *d {
            continue;
        }
        let (i, j) = (min(i, j), max(i, j));
        if let (Some((pi, si, di)), Some((pj, sj, dj))) = (parts[i].take(), parts[j].take()) {
            // Merge phylogeny i and j, adding theirs sizes. Store the result in the lower of the two indices.
            let inner_node_name = format!("inner_{}", inner_nodes_count);
            parts[i] = Some((
                Phylogeny::join(pi, *d / 2. - di, pj, *d / 2. - dj, &inner_node_name),
                si + sj,
                *d / 2.,
            ));
            inner_nodes_count += 1;

            // Update the distances to other nodes
            for k in 0..n {
                if k == i || k == j {
                    continue;
                }
                let dk =
                    (si as f64 * distances[i][k] + sj as f64 * distances[j][k]) / (si + sj) as f64;
                distances[i][k] = dk;
                distances[k][i] = dk;
                heap.push(Reverse((NotNan::new(dk).unwrap().into(), k, i)));
            }
        }
    }

    parts[0].take().unwrap().0
}
