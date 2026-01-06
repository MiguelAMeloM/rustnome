use crate::alignment::needleman_wunsch::align_profile;
use crate::alignment::Scoring;
use crate::sequences::profile::Profile;
use crate::sequences::DNuc;
use crate::upgma::tree::TreeNode;
use std::collections::HashMap;

fn create_profile<T: AsRef<[DNuc]>>(seq: T) -> Vec<Profile> {
    let mut prof = Vec::new();
    for nuc in seq.as_ref() {
        prof.push(nuc.prof());
    }
    prof
}

fn iter_step(
    nodes: &mut HashMap<usize, TreeNode>,
    matrix: &mut HashMap<usize, HashMap<usize, (f64, Vec<Profile>)>>,
    scoring: Scoring,
) {
    let (mut i, mut j) = (0, 0);
    let mut min_dist = f64::MAX;
    for (a, row) in matrix.iter() {
        for (b, (score, _)) in row {
            if *score < min_dist {
                min_dist = *score;
                i = *a;
                j = *b;
            }
        }
    }
    let node1 = nodes.remove(&i).unwrap();
    let node2 = nodes.remove(&j).unwrap();
    let (_, prof) = matrix.get(&i).unwrap().get(&j).unwrap();
    let new_node = node1.merge(node2, min_dist, prof.clone());
    matrix.remove(&i);
    matrix.remove(&j);
    for (_, row) in matrix.iter_mut() {
        row.remove(&j);
        row.remove(&i);
    }
    let max_idx = *nodes.keys().max().unwrap();
    let mut new_scores = HashMap::new();
    for (i, node) in nodes.iter() {
        let (score, new_prof) = align_profile(new_node.prof(), node.prof(), scoring);
        new_scores.insert(*i, (score, new_prof));
    }
    matrix.insert(max_idx + 1, new_scores);
    nodes.insert(max_idx + 1, new_node);
}

/// Constructs a UPGMA tree from a set of sequences using a given scoring scheme.
///
/// # Type Parameters
/// - `T`: A type that can be converted to a slice of `DNuc` (`AsRef<[DNuc]>`). This allows
///   the function to accept both vectors and slices of nucleotide sequences.
///
/// # Arguments
/// - `sequences`: A vector of tuples `(name, sequence)`, where `name` is a `String` identifying
///   the sequence, and `sequence` is of type `T` (convertible to a slice of `DNuc`).
/// - `scoring`: A `Scoring` object defining the parameters used to compute distances between sequences.
///
/// # Returns
/// - Returns a `TreeNode` representing the root of the constructed UPGMA tree.
///
/// # Behavior
/// - The function computes pairwise distances between sequences using the provided `scoring`
///   and progressively merges the closest pair of nodes until a single root node remains.
/// - The `TreeNode` structure contains either leaves (original sequences) or internal nodes
///   representing merged profiles and their associated distance.
///
/// # Example
/// ```rust
/// let sequences = vec![
///     ("seq1".to_string(), vec![DNuc::A, DNuc::C, DNuc::G]),
///     ("seq2".to_string(), vec![DNuc::A, DNuc::G, DNuc::G]),
/// ];
/// let scoring = Scoring::default();
/// let tree = create_tree(sequences, scoring);
/// ```
pub fn create_tree<T: AsRef<[DNuc]>>(sequences: Vec<(String, T)>, scoring: Scoring) -> TreeNode {
    let mut nodes = HashMap::new();
    let mut matrix = HashMap::new();
    for i in 0..sequences.len() {
        matrix.insert(i, HashMap::new());
        for j in i + 1..sequences.len() {
            let (_, node1) = sequences.get(i).unwrap();
            let (_, node2) = sequences.get(j).unwrap();
            matrix.get_mut(&i).unwrap().insert(
                j,
                align_profile(create_profile(node1), create_profile(node2), scoring),
            );
        }
    }
    for (i, (name, seq)) in sequences.into_iter().enumerate() {
        let prof = create_profile(seq);
        nodes.insert(i, TreeNode::Leaf(prof, name));
    }
    while nodes.len() > 1 {
        iter_step(&mut nodes, &mut matrix, scoring);
    }
    let (_, tree) = nodes.drain().next().unwrap();
    tree
}
