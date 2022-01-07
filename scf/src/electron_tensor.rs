use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use std::io::{stdout, Write};
use std::ops::Index;

pub struct ElectronTensor {
    n_basis: usize,
    data: Vec<f64>,
}

impl ElectronTensor {
    pub fn from_basis(basis: &[ContractedGaussian]) -> Self {
        let n_basis = basis.len();
        let data = (0..n_basis)
            .flat_map(move |w| {
                (0..n_basis).flat_map(move |z| {
                    (0..n_basis).flat_map(move |y| {
                        if y % 5 == 0 {
                            let progress = w * n_basis.pow(3) + z * n_basis.pow(2) + y * n_basis;
                            let max = n_basis * (n_basis.pow(3) + n_basis.pow(2) + n_basis + 1);
                            print!(
                                "\rERI-Formation progess: {}/{} ({:0.3}%)",
                                progress,
                                max,
                                progress as f32 / max as f32 * 100.0
                            );
                            stdout().flush().expect("failed to flush console");
                        }
                        (0..n_basis).map(move |x| {
                            if x >= y && z >= w && x * (x + 1) / 2 + y >= z * (z + 1) / 2 + w {
                                // TODO: implement screening routines
                                ContractedGaussian::electron_repulsion_int(
                                    &basis[x], &basis[y], &basis[z], &basis[w],
                                )
                            } else {
                                f64::NAN
                            }
                        })
                    })
                })
            })
            .collect::<Vec<_>>();

        Self { n_basis, data }
    }
}

impl Index<(usize, usize, usize, usize)> for ElectronTensor {
    type Output = f64;

    fn index(&self, (x, y, z, w): (usize, usize, usize, usize)) -> &Self::Output {
        let (x, y) = if x > y { (x, y) } else { (y, x) };
        let (z, w) = if z > w { (z, w) } else { (w, z) };
        let (x, y, z, w) = if x * (x + 1) / 2 + y >= z * (z + 1) / 2 + w {
            (x, y, z, w)
        } else {
            (z, w, x, y)
        };
        &self.data[x + y * self.n_basis + z * self.n_basis.pow(2) + w * self.n_basis.pow(3)]
    }
}
