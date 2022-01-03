use std::fmt::{Display, Formatter};
use std::ptr::write;

#[derive(Copy, Clone, Debug)]
pub struct GaussianPrimitive {
    angular: [i32; 3],
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

impl Display for GaussianPrimitive {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:0.4}*", self.coefficient)?;
        for (i, angular) in self.angular.into_iter().enumerate() {
            if angular > 1 {
                write!(f, "{}^{}*", ['x', 'y', 'z'][i], angular)?;
            } else if angular > 0 {
                write!(f, "{}*", ['x', 'y', 'z'][i])?;
            }
        }
        write!(f, "exp(-{:0.4}*r^2)", self.exponent)
    }
}

#[test]
pub fn test_display() {
    let primitive = GaussianPrimitive::new([2, 1, 0], 0.1266712, 1.0);
    println!("{}", primitive);
}
