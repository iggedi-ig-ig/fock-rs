use nalgebra::Vector3;

pub struct GaussianPrimitive {
    angular: Vector3<i32>,
    exponent: f64,
    coefficient: f64,
}

impl GaussianPrimitive {
    /// see https://www.wikiwand.com/en/Gaussian_orbital
    fn norm([i, j, k]: [i32; 3], alpha: f64) -> f64 {
        (std::f64::consts::FRAC_2_PI * alpha).powi(3).sqrt().sqrt()
            * f64::sqrt(
                (8.0 * alpha).powi(i + j + k)
                    / ((i + 1..=2 * i).product::<i32>()
                        * (j + 1..=2 * j).product::<i32>()
                        * (k + 1..=2 * k).product::<i32>()) as f64,
            )
    }

    pub fn new(angular: Vector3<i32>, exponent: f64, coefficient: f64) -> Self {
        let norm = Self::norm(angular.into(), exponent);
        GaussianPrimitive {
            angular,
            exponent,
            coefficient: coefficient * norm,
        }
    }

    pub fn new_spherical(exponent: f64, coefficient: f64) -> Self {
        Self::new(Vector3::zeros(), exponent, coefficient)
    }

    pub fn angular(&self) -> Vector3<i32> {
        self.angular
    }
    pub fn exponent(&self) -> f64 {
        self.exponent
    }
    pub fn coefficient(&self) -> f64 {
        self.coefficient
    }
}
