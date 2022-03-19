use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use log::info;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Index;

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
        let mut data = HashMap::new();

        for w in 0..n_basis {
            let perc = w as f32 / n_basis as f32 * 100.0;
            info!("Integral quadruplet {w}/{n_basis} ({perc:.1}% done)",);
            for z in 0..n_basis {
                for y in 0..n_basis {
                    for x in 0..n_basis {
                        let index = IntegralIndex::new((x, y, z, w));
                        data.entry(index).or_insert_with(|| {
                            ContractedGaussian::electron_repulsion_int(
                                &basis[x], &basis[y], &basis[z], &basis[w],
                            )
                        });
                    }
                }
            }
        }

        Self { data }
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
