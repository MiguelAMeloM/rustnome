/// Scoring parameters used for sequence alignment.
///
/// This struct groups together the numerical values that define how an
/// alignment is scored. It is intended to be passed as a single parameter
/// to alignment algorithms, ensuring consistency and clarity.
///
/// The scoring scheme consists of:
///
/// - A *gap penalty*, applied when a nucleotide is aligned with a gap.
/// - A *match score*, applied when two identical nucleotides are aligned.
/// - A *mismatch penalty*, applied when two different nucleotides are aligned.
///
/// Higher scores indicate better alignments. Gap and mismatch values are
/// typically negative, while match values are usually positive.
///
/// The struct is `Copy` and `Clone`, allowing it to be cheaply reused during
/// dynamic programming and backtracking.
#[derive(Copy, Clone)]
pub struct Scoring {
    /// Penalty applied when introducing a gap in the alignment.
    pub gap: f64,
    /// Penalty applied when aligning two different nucleotides.
    pub nuc_mismatch: f64,
    /// Score applied when aligning two identical nucleotides.
    pub nuc_match: f64,
}

pub trait Score {
    fn score(&self, other: &Self, scoring: Scoring) -> f64;
}

/// Type alias for a pair of aligned sequences.
///
/// Each sequence is represented as a vector of `Option<T>`, where `Some(T)`
/// indicates an aligned element and `None` represents a gap.
pub type Alignment<T> = [Vec<Option<T>>; 2];

mod matrix {
    pub type Matrix = Vec<Vec<f64>>;
    pub fn new_matrix(m: usize, n: usize) -> Matrix {
        let matrix = vec![vec![0.0; n + 1]; m + 1];
        matrix
    }
}

/// Module implementing the Smith-Waterman local alignment algorithm.
pub mod smith_waterman {
    use crate::alignment::{Alignment, Score, Scoring};
    use crate::alignment::matrix::{new_matrix, Matrix};

    fn find_max(matrix: &Matrix) -> (usize, usize) {
        let (mut x, mut y) = (0, 0);
        let mut max = 0.0;
        for i in 0..matrix.len() {
            for j in 0..matrix[i].len() {
                if matrix[i][j] > max {
                    max = matrix[i][j];
                    x = i;
                    y = j;
                };
            };
        };
        (x, y)
    }

    fn fill_matrix<T: Score>(
        seq1: &Vec<T>,
        seq2: &Vec<T>,
        matrix: &mut Matrix,
        scoring: Scoring,
    ) {
        let m = seq1.len();
        let n = seq2.len();
        for i in 1..=m {
            for j in 1..=n {
                let hor = matrix[i - 1][j] + scoring.gap;
                let ver = matrix[i][j - 1] + scoring.gap;
                let diag= matrix[i - 1][j - 1] + seq1[i - 1].score(&seq2[j - 1], scoring);
                matrix[i][j] = hor.max(ver).max(diag).max(0.0);
            };
        };
    }

    fn local_backtrack<T: Clone + Score>(
        als: &mut Vec<Alignment<T>>,
        matrix: &Matrix,
        (mut al1, mut al2): (Vec<Option<T>>, Vec<Option<T>>),
        (i, j): (usize, usize),
        (seq1, seq2): (&Vec<T>, &Vec<T>),
        punctuation: Scoring,
    ) {
        if matrix[i][j] == 0.0 {
            al1.reverse();
            al2.reverse();
            als.push([al1, al2]);
            return;
        }
        let score = matrix[i][j];
        if i > 0 && score == matrix[i - 1][j] + punctuation.gap {
            al1.push(Some(seq1[i - 1].clone()));
            al2.push(None);
            local_backtrack(als, matrix, (al1.clone(), al2.clone()), (i - 1, j), (seq1, seq2), punctuation);
        };
        if j > 0 && score == matrix[i][j - 1] + punctuation.gap {
            al2.push(Some(seq2[j - 1].clone()));
            al1.push(None);
            local_backtrack(als, matrix, (al1.clone(), al2.clone()), (i, j - 1), (seq1, seq2), punctuation);
        };
        if (i > 0 && j > 0) && (matrix[i][j] == matrix[i-1][j-1] + seq1[i - 1].score(&seq2[j - 1], punctuation))
        {
            al1.push(Some(seq1[i - 1].clone()));
            al2.push(Some(seq2[j - 1].clone()));
            local_backtrack(als, matrix, (al1.clone(), al2.clone()), (i - 1, j - 1), (seq1, seq2), punctuation);
        };
    }

    /// Performs a local alignment of two sequences using Smith-Waterman.
    ///
    /// # Parameters
    ///
    /// - `target`: First sequence.
    /// - `source`: Second sequence.
    /// - `punctuation`: Scoring parameters.
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - The maximum alignment score.
    /// - A vector of all optimal local alignments.
    pub fn align<T, U>(target: U, source: U, punctuation: Scoring) -> (f64, Vec<Alignment<T>>)
    where
        T: Score + Clone,
        U: AsRef<Vec<T>>,
    {
        let tar = target.as_ref();
        let src = source.as_ref();
        let mut matrix = new_matrix(tar.len(), src.len());
        fill_matrix(&tar, &src, &mut matrix, punctuation);
        let al1 = Vec::new();
        let al2 = Vec::new();
        let (i, j) = find_max(&matrix);
        let score = matrix[i][j];
        let mut als = Vec::new();
        local_backtrack(&mut als, &matrix, (al1, al2), (i, j), (tar.as_ref(), src.as_ref()), punctuation);
        (score, als)
    }
}

/// Module implementing the Needleman-Wunsch global alignment algorithm.
pub mod needleman_wunsch {
    use crate::alignment::{Alignment, Score, Scoring};
    use crate::sequences::profile::Profile;
    use super::matrix::*;

    fn init_matrix(matrix: &mut Matrix, gap: f64) {
        let m = matrix.len();
        let n = matrix[0].len();
        for i in 0..m {
            matrix[i][0] = i as f64 * gap;
        }
        for j in 0..n {
            matrix[0][j] = j as f64 * gap;
        }
    }

    fn fill_matrix<T: Score>(
        seq1: &Vec<T>,
        seq2: &Vec<T>,
        matrix: &mut Matrix,
        scoring: Scoring,
    ) {
        let m = seq1.len();
        let n = seq2.len();
        for i in 1..=m {
            for j in 1..=n {
                let hor = matrix[i - 1][j] + scoring.gap;
                let ver = matrix[i][j - 1] + scoring.gap;
                let diag= matrix[i - 1][j - 1] + seq1[i - 1].score(&seq2[j - 1], scoring);
                matrix[i][j] = hor.max(ver).max(diag);
            };
        };
    }

    fn backtrack<T: Clone + Score>(
        als: &mut Vec<Alignment<T>>,
        matrix: &Matrix,
        (mut al1, mut al2): (Vec<Option<T>>, Vec<Option<T>>),
        (i, j): (usize, usize),
        (seq1, seq2): (&Vec<T>, &Vec<T>),
        punctuation: Scoring,
    ) {
        if i == 0 && j == 0 {
            al1.reverse();
            al2.reverse();
            als.push([al1, al2]);
            return;
        }
        let score = matrix[i][j];
        if i > 0 && score == matrix[i - 1][j] + punctuation.gap {
            al1.push(Some(seq1[i - 1].clone()));
            al2.push(None);
            backtrack(als, matrix, (al1.clone(), al2.clone()), (i - 1, j), (seq1, seq2), punctuation);
        };
        if j > 0 && score == matrix[i][j - 1] + punctuation.gap {
            al2.push(Some(seq2[j - 1].clone()));
            al1.push(None);
            backtrack(als, matrix, (al1.clone(), al2.clone()), (i, j - 1), (seq1, seq2), punctuation);
        };
        if (i > 0 && j > 0) && (matrix[i][j] == matrix[i-1][j-1] + seq1[i - 1].score(&seq2[j - 1], punctuation))
        {
            al1.push(Some(seq1[i - 1].clone()));
            al2.push(Some(seq2[j - 1].clone()));
            backtrack(als, matrix, (al1.clone(), al2.clone()), (i - 1, j - 1), (seq1, seq2), punctuation);
        };
    }

    /// Performs a global alignment of two sequences using Needleman-Wunsch.
    ///
    /// # Parameters
    ///
    /// - `seq1`: First sequence.
    /// - `seq2`: Second sequence.
    /// - `punctuation`: Scoring parameters.
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - The optimal global alignment score.
    /// - A vector of all optimal global alignments.
    ///
    /// Gaps are represented as `None` and aligned elements as `Some(T)`.
    pub fn align<T, U>(seq1: U, seq2: U, scoring: Scoring) -> (f64, Vec<Alignment<T>>)
    where
        T: Score + Clone,
        U: AsRef<Vec<T>>,
    {
        let seq1 = seq1.as_ref();
        let seq2 = seq2.as_ref();
        let mut matrix = new_matrix(seq1.len(), seq2.len());
        init_matrix(&mut matrix, scoring.gap);
        fill_matrix(&seq1, &seq2, &mut matrix, scoring);
        let al1 = Vec::new();
        let al2 = Vec::new();
        let i = seq1.len();
        let j = seq2.len();
        let score = matrix[i][j];
        let mut als = Vec::new();
        backtrack(&mut als, &matrix, (al1, al2), (i, j), (seq1.as_ref(), seq2.as_ref()), scoring);
        (score, als)
    }


    /// Performs a global alignment of two profiles using Needleman-Wunsch.
    ///
    /// # Parameters
    ///
    /// - `prof1`: First sequence.
    /// - `prof2`: Second sequence.
    /// - `punctuation`: Scoring parameters.
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// - The optimal global alignment score.
    /// - The new profile representing the alignment.
    pub fn align_profile<T: AsRef<Vec<Profile>>>(prof1: T, prof2: T, scoring: Scoring) -> (f64, Vec<Profile>) {
        let prof1 = prof1.as_ref();
        let prof2 = prof2.as_ref();
        let mut matrix = new_matrix(prof1.len(), prof2.len());
        init_matrix(&mut matrix, scoring.gap);
        fill_matrix(&prof1, &prof2, &mut matrix, scoring);
        let mut result = Vec::new();
        let (mut i, mut j) = (prof1.len(), prof2.len());
        let score = matrix[i][j];
        while i > 0 || j > 0 {
            let score = matrix[i][j];
            if score == matrix[i - 1][j] + scoring.gap {
                result.push(prof1[i - 1].clone());
                i -= 1;
            } else if score == matrix[i][j - 1] + scoring.gap {
                result.push(prof2[j - 1].clone());
                j -= 1;
            } else {
                result.push(prof1[i-1].add(&prof2[j-1]))
            }
        };
        result.reverse();
        (score, result)
    }
}
