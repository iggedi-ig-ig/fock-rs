use std::ops::Index;

#[derive(Debug)]
pub struct ElectronRepulsionTensor {
    size: usize,
    data: Vec<f64>,
}

// TODO: use symmetry for better performance
impl ElectronRepulsionTensor {
    pub fn new_zeros(size: usize) -> Self {
        Self {
            size,
            data: vec![0.0; size * size * size * size],
        }
    }

    pub fn from_fn(size: usize, mut func: impl FnMut(usize, usize, usize, usize) -> f64) -> Self {
        let total = size * size * size * size;
        Self {
            size,
            data: (0..total)
                .map(|index| {
                    let x = index % size;
                    let y = ((index - x) / size) % size;
                    let z = ((index - y * size - x) / size.pow(2)) % size;
                    let w = (index - z * size.pow(2) - y * size - x) / size.pow(3);

                    func(x, y, z, w)
                })
                .collect::<Vec<_>>(),
        }
    }
}

impl Index<(usize, usize, usize, usize)> for ElectronRepulsionTensor {
    type Output = f64;

    fn index(&self, (x, y, z, w): (usize, usize, usize, usize)) -> &Self::Output {
        let index = x + y * self.size + z * self.size.pow(2) + w * self.size.pow(3);
        &self.data[index]
    }
}
