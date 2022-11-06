use nalgebra::Vector3;

#[derive(Copy, Clone, Debug)]
pub struct GaussianPrimitive {
    pub(crate) angular: [i32; 3],
    pub(crate) exponent: f64,
    pub coefficient: f64,
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

    pub fn new_unnormalized(angular: [i32; 3], exponent: f64, coefficient: f64) -> Self {
        GaussianPrimitive {
            angular,
            exponent,
            coefficient,
        }
    }

    pub fn new(angular: [i32; 3], exponent: f64, coefficient: f64) -> Self {
        let norm = Self::norm(angular, exponent);
        GaussianPrimitive {
            angular,
            exponent,
            coefficient: coefficient * norm,
        }
    }

    pub fn new_spherical(exponent: f64, coefficient: f64) -> Self {
        Self::new([0; 3], exponent, coefficient)
    }

    #[inline(always)]
    pub fn product_center(
        _a @ &Self {
            exponent: a_exp, ..
        }: &Self,
        a_pos: &Vector3<f64>,
        _b @ &Self {
            exponent: b_exp, ..
        }: &Self,
        b_pos: &Vector3<f64>,
    ) -> Vector3<f64> {
        (a_exp * a_pos + b_exp * b_pos) / (a_exp + b_exp)
    }

    pub fn angular(&self) -> [i32; 3] {
        self.angular
    }
    pub fn exponent(&self) -> f64 {
        self.exponent
    }
    pub fn coefficient(&self) -> f64 {
        self.coefficient
    }
}
