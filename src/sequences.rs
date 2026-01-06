use std::vec::IntoIter;
use crate::alignment::{Score, Scoring};
use crate::sequences::profile::Profile;

/// A DNA nucleotide.
///
/// This enum represents one of the four canonical nucleotides found in DNA.
/// Each variant corresponds to a specific nitrogenous base:
///
/// - `A`: Adenine
/// - `C`: Cytosine
/// - `T`: Thymine
/// - `G`: Guanine
///
/// `DNuc` is designed to be lightweight and copyable, making it suitable
/// for use in large DNA sequences without additional allocation costs.
#[derive(Copy, Clone)]
pub enum DNuc {
    /// Adenine, a purine base that pairs with Thymine.
    A,
    /// Cytosine, a pyrimidine base that pairs with Guanine.
    C,
    /// Thymine, a pyrimidine base that pairs with Adenine.
    T,
    /// Guanine, a purine base that pairs with Cytosine.
    G
}

impl DNuc {
    pub fn prof(&self) -> Profile {
        match self {
            DNuc::A => Profile::new(1.0, 0.0, 0.0, 0.0),
            DNuc::C => Profile::new(0.0, 1.0, 0.0, 0.0),
            DNuc::T => Profile::new(0.0, 0.0, 1.0, 0.0),
            DNuc::G => Profile::new(0.0, 0.0, 0.0, 1.0)
        }
    }
}

impl From<char> for DNuc {
    fn from(c: char) -> DNuc {
        match c {
            'a' | 'A' => DNuc::A,
            'c' | 'C' => DNuc::C,
            't' | 'T' => DNuc::T,
            'g' | 'G' => DNuc::G,
            _ => unreachable!()
        }
    }
}

impl PartialEq for DNuc {
    fn eq(&self, other: &DNuc) -> bool {
        match self {
            &DNuc::A => {
                if let &DNuc::A = other {
                    true
                } else {
                    false
                }
            },
            &DNuc::C => {
                if let &DNuc::C = other {
                    true
                } else {
                    false
                }
            },
            &DNuc::T => {
                if let &DNuc::T = other {
                    true
                } else {
                    false
                }
            },
            &DNuc::G => {
                if let &DNuc::G = other {
                    true
                } else {
                    false
                }
            }
        }
    }
}

impl Score for DNuc {
    fn score(&self, other: &DNuc, scoring: Scoring) -> f64 {
        if self == other {
            scoring.nuc_match
        } else {
            scoring.nuc_mismatch
        }
    }
}

/// A DNA sequence.
///
/// This struct represents a DNA molecule as an ordered sequence of
/// nucleotides (`DNuc`). Internally, the sequence is stored as a `Vec<DNuc>`,
/// preserving the original order of the bases.
///
/// The `DNA` type provides safe and explicit conversions from string-based
/// representations, allowing biological sequences to be constructed from
/// textual data such as FASTA records or user input.
pub struct DNA(Vec<DNuc>);

impl AsRef<Vec<DNuc>> for DNA {
    fn as_ref(&self) -> &Vec<DNuc> {
        &self.0
    }
}

impl From<&str> for DNA {
    fn from(s: &str) -> DNA {
        let mut seq = Vec::new();
        for c in s.chars() {
            seq.push(DNuc::from(c));
        };
        DNA(seq)
    }
}

impl From<String> for DNA {
    fn from(s: String) -> DNA {
        let mut seq = Vec::new();
        for c in s.chars() {
            seq.push(DNuc::from(c));
        };
        DNA(seq)
    }
}

impl IntoIterator for DNA {
    type Item = DNuc;
    type IntoIter = IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}


pub mod profile {
    use crate::alignment::{Score, Scoring};

    #[derive(Copy, Clone)]
    pub struct Profile(f64, f64, f64, f64);
    
    
    impl Profile {
        
        pub fn new(a: f64, c: f64, t: f64, g: f64) -> Profile {
            Profile(a, c, g, t)
        }
        fn norm(&self) -> f64 {
            (self.0.powi(2) + self.1.powi(2) + self.2.powi(2) + self.3.powi(2)).sqrt()
        }
        
        fn dot(&self, other: &Self) -> f64 {
            self.0*other.0 + self.1*other.1 + self.2*other.2 + self.3*other.3
        }
        
        fn similarity(&self, other: &Self) -> f64 {
            self.dot(other) / self.norm() / other.norm()
        }
        
        pub fn add(&self, other: &Self) -> Self {
            Profile(
                (self.0 + other.0)/2.0,
                (self.1 + other.1)/2.0,
                (self.2 + other.2)/2.0,
                (self.3 + other.3)/2.0
            )
        }
    }
    
    impl Score for Profile {
        fn score(&self, other: &Self, scoring: Scoring) -> f64 {
            let s = self.similarity(other);
            scoring.nuc_mismatch + s*(scoring.nuc_match - scoring.nuc_mismatch)
        }
    }
}

