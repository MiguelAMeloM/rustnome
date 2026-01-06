use std::collections::HashMap;

const LOCI_NUMBER: f64 = 3.0;
pub type Alleles = [usize; 2];

#[derive(Eq, PartialEq, Hash)]
pub enum Locus {
    AA1,
    AA2,
    X34,
}

pub type Genotipe = HashMap<Locus, Alleles>;

enum Coincidence {
    Cero,
    One,
    Two,
}

fn compare_aleles(a: &Alleles, b: &Alleles) -> Coincidence {
    if a[0] == b[0] && a[1] == b[1] {
        Coincidence::Two
    } else if a[0] == b[0] || a[1] == b[1] {
        Coincidence::One
    } else {
        Coincidence::Cero
    }
}

pub type InheritancePattern = (f64, f64, f64);

pub fn compare_genotypes(a: &Genotipe, b: &Genotipe) -> InheritancePattern {
    let (mut cero, mut one, mut two) = (0.0, 0.0, 0.0);
    for locus in a.keys() {
        let g1 = a.get(locus).unwrap();
        let g2 = b.get(locus).unwrap();
        match compare_aleles(g1, g2) {
            Coincidence::Cero => cero += 1.0,
            Coincidence::One => one += 1.0,
            Coincidence::Two => two += 1.0,
        }
    }
    (cero / LOCI_NUMBER, one / LOCI_NUMBER, two / LOCI_NUMBER)
}
