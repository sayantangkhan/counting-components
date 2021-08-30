//! Python module written in Rust to resolve via surgery an intersection of curves

#![feature(map_first_last)]

use gcd::Gcd;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::PyObjectProtocol;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashSet};

create_exception!(counting_components, PermutationException, PyException);

/// Enum describing possible errors when creating a signed permutation or multiple strands
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum PermutationError {
    /// Not a valid permutation on {0,1,...,n-1}
    InvalidPermutation,
    /// Flipset not a subset of {0,1,...,n-1}
    InvalidFlipset,
    /// Strand must be Transverse ('t') or a PermutationDirection ('p')
    InvalidStrandType,
}

impl std::convert::From<PermutationError> for PyErr {
    fn from(err: PermutationError) -> PyErr {
        match err {
            PermutationError::InvalidPermutation => {
                PermutationException::new_err("Invalid permutation")
            }
            PermutationError::InvalidFlipset => PermutationException::new_err("Invalid flip set"),
            PermutationError::InvalidStrandType => {
                PermutationException::new_err("Invalid strand type: only 't' and 'p' allowed")
            }
        }
    }
}

/// Permutation and flip data
#[pyclass]
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct SignedPermutation {
    permutation: Vec<usize>,
    flip_set: HashSet<usize>,
}

#[pymethods]
impl SignedPermutation {
    #[new]
    #[args(flips = "vec![]")]
    fn new(permutation: Vec<usize>, flips: Vec<usize>) -> PyResult<Self> {
        let length = permutation.len();
        let mut perm_vector = vec![length; length];
        let mut flip_set = HashSet::new();

        for (index, value) in permutation.into_iter().enumerate() {
            if value >= length {
                return Err(PermutationError::InvalidPermutation.into());
            }
            if perm_vector[value] != length {
                return Err(PermutationError::InvalidPermutation.into());
            }
            perm_vector[value] = index;
        }

        for value in flips.into_iter() {
            if value >= length {
                return Err(PermutationError::InvalidFlipset.into());
            }
            flip_set.insert(value);
        }

        Ok(Self {
            permutation: perm_vector,
            flip_set,
        })
    }

    #[call]
    fn __call__(&self, input: usize) -> PyResult<(usize, usize)> {
        if input >= self.permutation.len() {
            return Err(PermutationError::InvalidPermutation.into());
        }
        if self.flip_set.contains(&input) {
            Ok((self.permutation[input], 1))
        } else {
            Ok((self.permutation[input], 0))
        }
    }
}

#[pyproto]
impl PyObjectProtocol for SignedPermutation {
    fn __repr__(&self) -> PyResult<String> {
        let mut s = String::new();
        s.push('[');
        for (input, output) in self.permutation.iter().enumerate() {
            if self.flip_set.contains(&input) {
                s.push_str(&format!("{0} -> -{1}", input, output));
            } else {
                s.push_str(&format!("{0} -> {1}", input, output));
            }
            if input < self.permutation.len() - 1 {
                s.push_str(", ");
            }
        }
        s.push(']');
        Ok(s)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, PartialOrd, Ord)]
enum Strand {
    Transverse(usize),
    PermutationDirection(usize, usize),
}

/// Python class to represent a strand
#[pyclass]
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, PartialOrd, Ord)]
struct PyStrand {
    strand: Strand,
}

#[pymethods]
impl PyStrand {
    #[new]
    #[args(n = "0")]
    fn new(type_of_strand: char, m: usize, n: usize) -> PyResult<Self> {
        match type_of_strand {
            't' => Ok(Self {
                strand: Strand::Transverse(m),
            }),
            'p' => Ok(Self {
                strand: Strand::PermutationDirection(m, n),
            }),
            _ => Err(PermutationError::InvalidStrandType.into()),
        }
    }
}

#[pyproto]
impl PyObjectProtocol for PyStrand {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("{:?}", self.strand))
    }
}

/// Returns next major strand and info about whether it flipped
#[pyfunction]
fn get_next_major_strand(
    perm: &SignedPermutation,
    m: usize,
    n: usize,
    strand: PyStrand,
) -> (PyStrand, usize) {
    let mut flipped = 0;
    let out_strand: PyStrand = match strand.strand {
        Strand::PermutationDirection(mut perm_index, mut copy_index) => {
            if perm.flip_set.contains(&perm_index) {
                copy_index = m - copy_index - 1;
                flipped = 1;
            }

            perm_index = perm.permutation[perm_index];
            let mut absolute_index = m * perm_index + copy_index;

            if absolute_index + n < perm.permutation.len() * m {
                absolute_index += n;
                perm_index = absolute_index / m;
                copy_index = absolute_index % m;
                PyStrand {
                    strand: Strand::PermutationDirection(perm_index, copy_index),
                }
            } else {
                let transverse_index = perm.permutation.len() * m - absolute_index - 1;
                PyStrand {
                    strand: Strand::Transverse(transverse_index),
                }
            }
        }
        Strand::Transverse(index) => {
            if index + perm.permutation.len() * m < n {
                PyStrand {
                    strand: Strand::Transverse(index + perm.permutation.len() * m),
                }
            } else {
                let absolute_index = n - index - 1;
                let perm_index = absolute_index / m;
                let copy_index = absolute_index % m;
                PyStrand {
                    strand: Strand::PermutationDirection(perm_index, copy_index),
                }
            }
        }
    };
    (out_strand, flipped)
}

/// Determines if a given (perm, m, n) triple has only one component and outputs the orientability
/// Can I make this algorithm parallel?
#[pyfunction]
fn has_one_component(perm: &SignedPermutation, m: usize, n: usize) -> (bool, usize) {
    let expected_orbit_length = m * perm.permutation.len() + n;
    let mut actual_orbit_length = 1;

    let starting_strand = PyStrand {
        strand: Strand::Transverse(0),
    };
    let mut next_strand_with_orientability = get_next_major_strand(perm, m, n, starting_strand);
    let mut next_strand = next_strand_with_orientability.0;
    let mut orientability = next_strand_with_orientability.1;
    while next_strand != starting_strand {
        next_strand_with_orientability = get_next_major_strand(perm, m, n, next_strand);
        next_strand = next_strand_with_orientability.0;
        orientability = (orientability + next_strand_with_orientability.1) % 2;
        actual_orbit_length += 1;
    }

    (expected_orbit_length == actual_orbit_length, orientability)
}

/// Count components with orientability: ouputs a tuple indicating the number of two-sided and one-sided components
/// Can I make this parallel?
#[pyfunction]
fn count_components_with_orientability(
    perm: &SignedPermutation,
    m: usize,
    n: usize,
) -> (usize, usize) {
    let mut two_sided_components = 0;
    let mut one_sided_components = 0;

    let mut strands = BTreeSet::new();
    for i in 0..n {
        strands.insert(PyStrand {
            strand: Strand::Transverse(i),
        });
    }
    for j in 0..perm.permutation.len() {
        for k in 0..m {
            strands.insert(PyStrand {
                strand: Strand::PermutationDirection(j, k),
            });
        }
    }

    while !strands.is_empty() {
        let first_strand = strands.pop_first().unwrap();
        let next_strand_with_orientability = get_next_major_strand(perm, m, n, first_strand);
        let mut orientability = next_strand_with_orientability.1;
        let mut next_strand = next_strand_with_orientability.0;
        while next_strand != first_strand {
            let next_strand_with_orientability = get_next_major_strand(perm, m, n, next_strand);
            strands.remove(&next_strand);
            orientability += next_strand_with_orientability.1;
            next_strand = next_strand_with_orientability.0;
        }
        if orientability % 2 == 0 {
            two_sided_components += 1;
        } else {
            one_sided_components += 1;
        }
    }

    (two_sided_components, one_sided_components)
}

/// Function to count components of all (m,n) pairs up to a complexity in parallel
#[pyfunction]
fn count_components_upto_complexity(
    perm: &SignedPermutation,
    complexity: usize,
) -> Vec<((usize, usize), (usize, usize))> {
    (2..complexity)
        .into_par_iter()
        .flat_map(|k| {
            (1..k)
                .into_par_iter()
                .filter(move |n| k.gcd_binary(*n) == 1)
                .map(move |n| {
                    let m = k - n;
                    ((m, n), count_components_with_orientability(perm, m, n))
                })
        })
        .collect()
}

/// Function to list only two-sided multicurves up to a given complexity
#[pyfunction]
fn two_sided_multicurves_upto_complexity(
    perm: &SignedPermutation,
    complexity: usize,
) -> Vec<(usize, usize)> {
    (2..complexity)
        .into_par_iter()
        .flat_map(|k| {
            (1..k)
                .into_par_iter()
                .filter(move |n| k.gcd_binary(*n) == 1)
                .map(move |n| {
                    let m = k - n;
                    ((m, n), count_components_with_orientability(perm, m, n))
                })
                .filter(|(_, (_, o))| *o == 0)
                .map(|(a, _)| a)
        })
        .collect()
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn counting_components(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<SignedPermutation>()?;
    m.add_class::<PyStrand>()?;
    m.add_function(wrap_pyfunction!(get_next_major_strand, m)?)?;
    m.add_function(wrap_pyfunction!(has_one_component, m)?)?;
    m.add_function(wrap_pyfunction!(count_components_with_orientability, m)?)?;
    m.add_function(wrap_pyfunction!(count_components_upto_complexity, m)?)?;
    m.add_function(wrap_pyfunction!(two_sided_multicurves_upto_complexity, m)?)?;
    m.add(
        "PermutationException",
        _py.get_type::<PermutationException>(),
    )?;

    Ok(())
}
