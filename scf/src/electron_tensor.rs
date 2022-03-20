use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use log::debug;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Index;
use std::sync::{Arc, Mutex};

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct IntegralIndex {
    x: usize,
    y: usize,
    z: usize,
    w: usize,
}

impl IntegralIndex {
    pub fn new(index: (usize, usize, usize, usize)) -> Self {
        let (x, y, z, w) = Self::correct_order(index);
        Self { x, y, z, w }
    }

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

pub struct ElectronTensor {
    data: HashMap<IntegralIndex, f64>,
}

impl ElectronTensor {
    pub fn from_basis(basis: &[ContractedGaussian]) -> Self {
        let n_basis = basis.len();
        let n_integrals = (n_basis.pow(4) + 8 + 1) / 8;
        let data = Arc::new(Mutex::new(HashMap::new()));

        (0..n_basis).into_par_iter().for_each(|w| {
            (w..n_basis).for_each(|z| {
                (0..n_basis).for_each(|y| {
                    (0..n_basis).for_each(|x| {
                        let index = IntegralIndex::new((x, y, z, w));

                        let lock = data.lock().unwrap();
                        if !lock.contains_key(&index) {
                            // drop lock so other threads can continue working
                            drop(lock);

                            let integral = ContractedGaussian::electron_repulsion_int(
                                &basis[x], &basis[y], &basis[z], &basis[w],
                            );

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
