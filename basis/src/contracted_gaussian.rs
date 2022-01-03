use crate::primitives::GaussianPrimitive;
use crate::utils::{coulomb_auxiliary, hermite_expansion};
use crate::BasisFunction;
use nalgebra::Vector3;

#[derive(Clone, Debug)]
pub struct ContractedGaussian {
    position: Vector3<f64>,
    primitives: Vec<GaussianPrimitive>,
}

impl ContractedGaussian {
    pub fn new(position: Vector3<f64>, primitives: Vec<GaussianPrimitive>) -> Self {
        ContractedGaussian {
            position,
            primitives,
        }
    }
    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn primitives(&self) -> &[GaussianPrimitive] {
        &self.primitives
    }
}

impl GaussianPrimitive {
    /// helper function to add angular momentum
    fn _add_ijk(primitive: &GaussianPrimitive, [i, j, k]: [i32; 3]) -> GaussianPrimitive {
        let [a, b, c] = primitive.angular();
        GaussianPrimitive::new_unnormalized(
            [a + i, b + j, c + k],
            primitive.exponent(),
            primitive.coefficient(),
        )
    }

    /// overlap integral between two primitive gaussians
    fn _overlap(a: &GaussianPrimitive, b: &GaussianPrimitive, diff: &Vector3<f64>) -> f64 {
        let [l1, m1, n1] = a.angular();
        let [l2, m2, n2] = b.angular();

        let s1 = hermite_expansion(l1, l2, 0, diff.x, a.exponent(), b.exponent());
        let s2 = hermite_expansion(m1, m2, 0, diff.y, a.exponent(), b.exponent());
        let s3 = hermite_expansion(n1, n2, 0, diff.z, a.exponent(), b.exponent());

        s1 * s2
            * s3
            * (std::f64::consts::PI / (a.exponent() + b.exponent()))
                .powi(3)
                .sqrt()
    }

    /// kinetic energy integral between two primitive gaussians
    fn _kinetic(a: &GaussianPrimitive, b: &GaussianPrimitive, diff: &Vector3<f64>) -> f64 {
        let [l2, m2, n2] = b.angular();

        let term0 = b.exponent() * (2 * (l2 + m2 + n2) + 3) as f64 * Self::_overlap(a, b, diff);
        let term1 = -2.0
            * b.exponent().powi(2)
            * (Self::_overlap(a, &Self::_add_ijk(b, [2, 0, 0]), diff)
                + Self::_overlap(a, &Self::_add_ijk(b, [0, 2, 0]), diff)
                + Self::_overlap(a, &Self::_add_ijk(b, [0, 0, 2]), diff));
        let term2 =
            -0.5 * (l2 * (l2 - 1)) as f64 * Self::_overlap(a, &Self::_add_ijk(a, [-2, 0, 0]), diff)
                + (m2 * (m2 - 1)) as f64 * Self::_overlap(a, &Self::_add_ijk(a, [0, -2, 0]), diff)
                + (n2 * (n2 - 1)) as f64 * Self::_overlap(a, &Self::_add_ijk(a, [0, 0, -2]), diff);
        term0 + term1 + term2
    }

    /// nuclear attraction integral between two primitive gaussians
    fn _nuclear_attraction(
        a: &GaussianPrimitive,
        b: &GaussianPrimitive,
        diff: &Vector3<f64>,
        prod_center: &Vector3<f64>,
        nucleus_pos: &Vector3<f64>,
        nucleus_charge: f64,
    ) -> f64 {
        let [l1, m1, n1] = a.angular();
        let [l2, m2, n2] = b.angular();

        let p = a.exponent() + b.exponent();

        let diff_comp_nucleus = prod_center - nucleus_pos;
        let dist_sq_comp_nucleus = diff_comp_nucleus.norm_squared();

        -nucleus_charge * std::f64::consts::TAU / p
            * (0..=l1 + l2)
                .flat_map(|t| {
                    (0..=m1 + m2).flat_map(move |u| {
                        (0..=n1 + n2).map(move |v| {
                            hermite_expansion(l1, l2, t, diff.x, a.exponent(), b.exponent())
                                * hermite_expansion(m1, m2, u, diff.y, a.exponent(), b.exponent())
                                * hermite_expansion(n1, n2, v, diff.z, a.exponent(), b.exponent())
                                * coulomb_auxiliary(
                                    t,
                                    u,
                                    v,
                                    0,
                                    p,
                                    &diff_comp_nucleus,
                                    dist_sq_comp_nucleus,
                                )
                        })
                    })
                })
                .sum::<f64>()
    }

    // TODO: check if this is correct
    #[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
    /// four-electron integrals between four primitive gaussians
    fn _electron_repulsion(
        a: &GaussianPrimitive,
        b: &GaussianPrimitive,
        c: &GaussianPrimitive,
        d: &GaussianPrimitive,
        diff_ab: &Vector3<f64>,
        diff_cd: &Vector3<f64>,
        comp_diff: &Vector3<f64>,
        comp_dist_sq: f64,
    ) -> f64 {
        let [l1, m1, n1]: [i32; 3] = a.angular();
        let [l2, m2, n2]: [i32; 3] = b.angular();
        let [l3, m3, n3]: [i32; 3] = c.angular();
        let [l4, m4, n4]: [i32; 3] = d.angular();

        let a = a.exponent();
        let b = b.exponent();
        let c = c.exponent();
        let d = d.exponent();

        let p = a + b;
        let q = c + d;

        let alpha = p * q / (p + q);
        2.0 * std::f64::consts::PI.powi(5).sqrt()
            * (p * q * (p + q).sqrt()).recip()
            * (0..=l1 + l2)
                .flat_map(move |t1| {
                    (0..=m1 + m2).flat_map(move |u1| {
                        (0..=n1 + n2).flat_map(move |v1| {
                            (0..=l3 + l4).flat_map(move |t2| {
                                (0..=m3 + m4).flat_map(move |u2| {
                                    (0..=n3 + n4).map(move |v2| {
                                        (-1.0f64).powi(t2 + u2 + v2)
                                            * hermite_expansion(l1, l2, t1, diff_ab.x, a, b)
                                            * hermite_expansion(m1, m2, u1, diff_ab.y, a, b)
                                            * hermite_expansion(n1, n2, v1, diff_ab.z, a, b)
                                            * hermite_expansion(l3, l4, t2, diff_cd.x, c, d)
                                            * hermite_expansion(m3, m4, u2, diff_cd.y, c, d)
                                            * hermite_expansion(n3, n4, v2, diff_cd.z, c, d)
                                            * coulomb_auxiliary(
                                                t1 + t2,
                                                u1 + u2,
                                                v1 + v2,
                                                0,
                                                alpha,
                                                &comp_diff,
                                                comp_dist_sq,
                                            )
                                    })
                                })
                            })
                        })
                    })
                })
                .sum::<f64>()
    }
}

impl BasisFunction for ContractedGaussian {
    fn overlap_int(a: &Self, b: &Self) -> f64 {
        let m = a.primitives.len();
        let n = b.primitives.len();
        let diff = b.position() - a.position();
        (0..m)
            .flat_map(move |i| {
                (0..n).map(move |j| {
                    a.primitives[i].coefficient()
                        * b.primitives[j].coefficient()
                        * GaussianPrimitive::_overlap(&a.primitives[i], &b.primitives[j], &diff)
                })
            })
            .sum::<f64>()
    }

    fn kinetic_int(a: &Self, b: &Self) -> f64 {
        let m = a.primitives.len();
        let n = b.primitives.len();
        let diff = b.position() - a.position();
        (0..m)
            .flat_map(move |i| {
                (0..n).map(move |j| {
                    a.primitives[i].coefficient()
                        * b.primitives[j].coefficient()
                        * GaussianPrimitive::_kinetic(&a.primitives[i], &b.primitives[j], &diff)
                })
            })
            .sum::<f64>()
    }

    fn nuclear_attraction_int(
        a: &Self,
        b: &Self,
        nucleus_pos: &Vector3<f64>,
        nucleus_charge: f64,
    ) -> f64 {
        let m = a.primitives.len();
        let n = b.primitives.len();
        let diff = b.position() - a.position();
        (0..m)
            .flat_map(move |i| {
                (0..n).map(move |j| {
                    let prod_center = a.primitives[i].exponent() * a.position()
                        + b.primitives[j].exponent() * b.position()
                            / (a.primitives[i].exponent() + b.primitives[j].exponent());

                    let a = &a.primitives[i];
                    let b = &b.primitives[j];

                    a.coefficient()
                        * b.coefficient()
                        * GaussianPrimitive::_nuclear_attraction(
                            a,
                            b,
                            &diff,
                            &prod_center,
                            nucleus_pos,
                            nucleus_charge,
                        )
                })
            })
            .sum::<f64>()
    }

    fn electron_repulsion_int(a: &Self, b: &Self, c: &Self, d: &Self) -> f64 {
        let m = a.primitives.len();
        let n = b.primitives.len();
        let o = c.primitives.len();
        let p = d.primitives.len();

        let diff_ab = b.position() - a.position();
        let diff_cd = d.position() - c.position();
        (0..m)
            .flat_map(move |i| {
                (0..n).flat_map(move |j| {
                    (0..o).flat_map(move |tau| {
                        (0..p).map(move |nu| {
                            let a_pos = a.position();
                            let a = &a.primitives[i];

                            let b_pos = b.position();
                            let b = &b.primitives[j];

                            let c_pos = c.position();
                            let c = &c.primitives[tau];

                            let d_pos = d.position();
                            let d = &d.primitives[nu];

                            let comp_ab = (a_pos * a.exponent() + b_pos * b.exponent())
                                / (a.exponent() + b.exponent());
                            let comp_cd = (c_pos * c.exponent() + d_pos * d.exponent())
                                / (c.exponent() + d.exponent());

                            let comp_diff = comp_cd - comp_ab;
                            let comp_dist_sq = comp_diff.norm_squared();

                            a.coefficient()
                                * b.coefficient()
                                * c.coefficient()
                                * d.coefficient()
                                * GaussianPrimitive::_electron_repulsion(
                                    a,
                                    b,
                                    c,
                                    d,
                                    &diff_ab,
                                    &diff_cd,
                                    &comp_diff,
                                    comp_dist_sq,
                                )
                        })
                    })
                })
            })
            .sum::<f64>()
    }

    fn evaluate(&self, at: &Vector3<f64>) -> f64 {
        let diff = self.position - at;
        let dist_sq = diff.norm_squared();
        self.primitives.iter().fold(0.0, |acc, gaussian| {
            let [i, j, k] = gaussian.angular();
            acc + gaussian.coefficient()
                * diff.x.powi(i)
                * diff.y.powi(j)
                * diff.z.powi(k)
                * f64::exp(-gaussian.exponent() * dist_sq)
        })
    }
}
