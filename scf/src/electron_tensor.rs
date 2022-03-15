use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use log::info;
use std::ops::Index;

#[derive(Debug)]
pub struct ElectronTensor {
    pub n_basis: usize,
    pub data: Vec<f64>,
}

impl ElectronTensor {
    pub fn from_basis(basis: &[ContractedGaussian]) -> Self {
        let n_basis = basis.len();
        let data = (0..n_basis)
            .flat_map(move |w| {
                info!(
                    "ERI-Formation progress: {w}/{n_basis} ({:0.3}%)",
                    w as f32 / n_basis as f32 * 100f32
                );
                (0..n_basis).flat_map(move |z| {
                    (0..n_basis).flat_map(move |y| {
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
