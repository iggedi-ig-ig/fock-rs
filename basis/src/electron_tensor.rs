use crate::contracted_gaussian::ContractedGaussian;
use crate::BasisFunction;
use log::debug;
use rayon::prelude::*;
use std::collections::hash_map::Entry::Vacant;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Index;
use std::sync::{Arc, Mutex};

/// An integral index used in the two-electron integrals of a basis set.
///
/// The index represents the four indices (x, y, z, w) used to calculate a two-electron integral:
///   int_{x,y,z,w} = int_{xy|zw} = <x y | z w>
///
/// Since two-electron integrals are symmetric in (xy) and (zw), this struct stores its indices
/// with the correct order: xy <= zw, reducing the total number of unique integrals to compute.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct IntegralIndex {
    x: usize,
    y: usize,
    z: usize,
    w: usize,
}

impl IntegralIndex {
    /// Creates a new integral index with the given indices.
    pub fn new(index: (usize, usize, usize, usize)) -> Self {
        let (x, y, z, w) = Self::correct_order(index);
        Self { x, y, z, w }
    }

    /// Returns the indices with the correct order, such that xy <= zw.
    pub fn correct_order(
        (x, y, z, w): (usize, usize, usize, usize),
    ) -> (usize, usize, usize, usize) {
        let (x, y) = if x > y { (x, y) } else { (y, x) };
        let (z, w) = if z > w { (z, w) } else { (w, z) };

        let xy = x * (x + 1) / 2 + y;
        let zw = z * (z + 1) / 2 + w;

        if xy > zw {
            (x, y, z, w)
        } else {
            (z, w, x, y)
        }
    }
}

/// An electron tensor representing electron-electron repulsion integrals between
/// four contracted Gaussian functions in a given basis set.
pub struct ElectronTensor {
    data: HashMap<IntegralIndex, f64>,
}

impl ElectronTensor {
    /// Constructs an `ElectronTensor` from the given basis set. Computes electron-electron
    /// repulsion integrals for each unique combination of four Gaussian functions in the basis
    /// set and stores them in a hashmap. This method utilizes parallel processing with the
    /// Rayon library to speed up computation time.
    ///
    /// # Arguments
    ///
    /// * `basis` - A slice of `ContractedGaussian` functions representing the basis set to use
    /// for computing electron-electron repulsion integrals.
    ///
    /// # Returns
    ///
    /// An `ElectronTensor` containing a hashmap of electron-electron repulsion integrals
    /// computed for each unique combination of four Gaussian functions in the given basis set.
    pub fn from_basis(basis: &[ContractedGaussian]) -> Self {
        // Initialize variables for computing the total number of integrals and a thread-safe
        // container for storing the resulting electron-electron repulsion integrals.
        let n_basis = basis.len();
        let n_integrals = (n_basis.pow(4) + 8 + 1) / 8;
        let mut diagonal = HashMap::with_capacity(n_integrals);

        // compute diagonal first
        (0..n_basis).for_each(|j| {
            (j..n_basis).for_each(|i| {
                let index = IntegralIndex::new((i, j, i, j));

                if let Vacant(e) = diagonal.entry(index) {
                    e.insert(ContractedGaussian::electron_repulsion_int(
                        &basis[i], &basis[j], &basis[i], &basis[j],
                    ));
                }
            })
        });

        let data = Arc::new(Mutex::new(diagonal.clone()));
        // Use parallel processing to compute electron-electron repulsion integrals for each
        // unique combination of four Gaussian functions in the basis set and store the result
        // in the hashmap.
        (0..n_basis).par_bridge().for_each(|w| {
            (w..n_basis).for_each(|z| {
                (0..n_basis).for_each(|y| {
                    (y..n_basis).for_each(|x| {
                        let index = IntegralIndex::new((x, y, z, w));

                        let mut lock = data.lock().unwrap();
                        if let Vacant(e) = lock.entry(index) {
                            // insert something into hashmap to prevent double computation
                            e.insert(f64::NAN);

                            // drop lock so other threads can continue working
                            drop(lock);

                            let diagonal_index_ij = IntegralIndex::new((x, y, x, y));
                            let diagonal_index_kl = IntegralIndex::new((z, w, z, w));

                            let estimate = f64::sqrt(
                                diagonal[&diagonal_index_ij] * diagonal[&diagonal_index_kl],
                            );

                            let integral = if estimate > 1e-6 {
                                ContractedGaussian::electron_repulsion_int(
                                    &basis[x], &basis[y], &basis[z], &basis[w],
                                )
                            } else {
                                0.0
                            };

                            data.lock().unwrap().insert(index, integral);
                        }
                    })
                })
            });
            let amount = data.lock().unwrap().len();
            debug!(
                "{amount} integrals computed, {:.1}% done",
                (amount as f32 / n_integrals as f32).min(1.0) * 100.0
            );
        });

        // Extract the hashmap from the thread-safe container and return an `ElectronTensor`.
        Self {
            data: Arc::try_unwrap(data).unwrap().into_inner().unwrap(),
        }
    }
}

impl Index<(usize, usize, usize, usize)> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        &self.data[&IntegralIndex::new(index)]
    }
}

impl Index<IntegralIndex> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: IntegralIndex) -> &Self::Output {
        &self.data[&index]
    }
}
