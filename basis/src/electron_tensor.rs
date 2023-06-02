use crate::contracted_gaussian::ContractedGaussian;
use crate::BasisFunction;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::hash::Hash;
use std::ops::Index;
use std::sync::RwLock;

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
        Self::new_unchecked(Self::correct_order(index))
    }

    /// Creates a new integral index, but doesn't correct their order
    pub fn new_unchecked((x, y, z, w): (usize, usize, usize, usize)) -> Self {
        Self { x, y, z, w }
    }

    /// Returns the indices with the correct order, such that xy <= zw.
    #[inline]
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

    pub fn linear(&self, size: usize) -> usize {
        self.w * size.pow(3) + self.z * size.pow(2) + self.y * size + self.x
    }
}

/// An electron tensor representing electron-electron repulsion integrals between
/// four contracted Gaussian functions in a given basis set.
pub struct ElectronTensor {
    data: Vec<f64>,
    /// side length
    size: usize,
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
        let data = (0..n_basis.pow(4))
            .map(|_| RwLock::new(None))
            .collect::<Vec<_>>();

        // compute diagonal first
        (0..n_basis).for_each(|j| {
            (j..n_basis).for_each(|i| {
                let index = IntegralIndex::new_unchecked((i, j, i, j));
                let linear = index.linear(n_basis);

                let _ = data[linear].write().unwrap().get_or_insert_with(|| {
                    ContractedGaussian::electron_repulsion_int(
                        &basis[i], &basis[j], &basis[i], &basis[j],
                    )
                });
            })
        });

        // Use parallel processing to compute electron-electron repulsion integrals for each
        // unique combination of four Gaussian functions in the basis set and store the result
        // in the hashmap.
        let pb = ProgressBar::new(((n_basis.pow(2) + 1) / 2).pow(2) as u64)
            .with_style(
                ProgressStyle::with_template("{prefix:>12.cyan.bold} [{bar:57}] {pos}/{len}")
                    .unwrap()
                    .progress_chars("=> "),
            )
            .with_prefix("compute ERIs");
        (0..n_basis).into_par_iter().for_each(|w| {
            (w..n_basis).for_each(|z| {
                (0..n_basis).for_each(|y| {
                    (y..n_basis).for_each(|x| {
                        let xy = x * (x + 1) / 2 + y;
                        let zw = z * (z + 1) / 2 + w;

                        // we know that x, y and z, w are always in the correct order (x <= y and z <= w)
                        // we thus only need to correct for hyper order
                        let index = IntegralIndex::new_unchecked(if xy > zw {
                            (x, y, z, w)
                        } else {
                            (z, w, x, y)
                        });
                        let linear = index.linear(n_basis);

                        let diagonal_index_ij =
                            IntegralIndex::new_unchecked((index.x, index.y, index.x, index.y));
                        let diagonal_index_kl =
                            IntegralIndex::new_unchecked((index.z, index.w, index.z, index.w));

                        let estimate = f64::sqrt(
                            data[diagonal_index_ij.linear(n_basis)]
                                .read()
                                .unwrap()
                                .expect("diagonal is set")
                                * data[diagonal_index_kl.linear(n_basis)]
                                    .read()
                                    .unwrap()
                                    .expect("diagonal is set"),
                        );

                        let lock = data[linear].read().unwrap();
                        if lock.is_none() {
                            drop(lock);

                            let integral = if estimate > 1e-6 {
                                ContractedGaussian::electron_repulsion_int(
                                    &basis[x], &basis[y], &basis[z], &basis[w],
                                )
                            } else {
                                0.0
                            };

                            *data[linear].write().unwrap() = Some(integral);
                        }
                    });

                    pb.inc((n_basis - y) as u64);
                })
            });
        });
        pb.finish();

        // Extract the hashmap from the thread-safe container and return an `ElectronTensor`.
        Self {
            data: data
                .into_iter()
                .map(|x| x.into_inner().unwrap())
                .map(|x| x.unwrap_or_default())
                .collect(),
            size: n_basis,
        }
    }
}

impl Index<(usize, usize, usize, usize)> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        let index = IntegralIndex::new(index);
        let linear = index.linear(self.size);
        &self.data[linear]
    }
}

impl Index<IntegralIndex> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: IntegralIndex) -> &Self::Output {
        let linear = index.linear(self.size);
        &self.data[linear]
    }
}
