use crate::primitives::GaussianPrimitive;
use crate::utils::{coulomb_auxiliary, hermite_expansion};
use crate::{BasisFunction, PointCharge};
use nalgebra::Vector3;
use smallvec::SmallVec;

#[derive(Clone, Debug)]
pub struct ContractedGaussian {
    position: Vector3<f64>,
    primitives: SmallVec<[GaussianPrimitive; Self::SMALL_VEC_SIZE]>,
}

impl ContractedGaussian {
    pub const SMALL_VEC_SIZE: usize = 6;

    pub fn new(
        position: Vector3<f64>,
        primitives: SmallVec<[GaussianPrimitive; Self::SMALL_VEC_SIZE]>,
    ) -> Self {
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
    fn _add_ijk(
        _primitive @ &GaussianPrimitive {
            angular: [a, b, c],
            exponent,
            coefficient,
        }: &GaussianPrimitive,
        [i, j, k]: [i32; 3],
    ) -> GaussianPrimitive {
        GaussianPrimitive::new_unnormalized([a + i, b + j, c + k], exponent, coefficient)
    }

    /// overlap integral between two primitive gaussians
    fn _overlap(
        _a @ &GaussianPrimitive {
            angular: [l1, m1, n1],
            exponent: a_exp,
            ..
        }: &GaussianPrimitive,
        _b @ &GaussianPrimitive {
            angular: [l2, m2, n2],
            exponent: b_exp,
            ..
        }: &GaussianPrimitive,
        diff: &Vector3<f64>,
    ) -> f64 {
        hermite_expansion([l1, l2, 0], diff.x, a_exp, b_exp)
            * hermite_expansion([m1, m2, 0], diff.y, a_exp, b_exp)
            * hermite_expansion([n1, n2, 0], diff.z, a_exp, b_exp)
            * (std::f64::consts::PI / (a_exp + b_exp)).powi(3).sqrt()
    }

    /// kinetic energy integral between two primitive gaussians
    fn _kinetic(
        a: &GaussianPrimitive,
        b @ &GaussianPrimitive {
            angular: [l2, m2, n2],
            exponent: b_exp,
            ..
        }: &GaussianPrimitive,
        diff: &Vector3<f64>,
    ) -> f64 {
        let overlap_step = |i, j, k| Self::_overlap(a, &Self::_add_ijk(b, [i, j, k]), diff);

        let term0 = b_exp * (2 * (l2 + m2 + n2) + 3) as f64 * Self::_overlap(a, b, diff);
        let term1 = -2.0
            * b.exponent().powi(2)
            * (overlap_step(2, 0, 0) + overlap_step(0, 2, 0) + overlap_step(0, 0, 2));
        let term2 = -0.5
            * ((l2 * (l2 - 1)) as f64 * overlap_step(-2, 0, 0)
                + (m2 * (m2 - 1)) as f64 * overlap_step(0, -2, 0)
                + (n2 * (n2 - 1)) as f64 * overlap_step(0, 0, -2));
        term0 + term1 + term2
    }

    /// nuclear attraction integral between two primitive gaussians
    pub fn _nuclear_attraction(
        _a @ &GaussianPrimitive {
            angular: [l1, m1, n1],
            exponent: a,
            ..
        }: &GaussianPrimitive,
        _b @ &GaussianPrimitive {
            angular: [l2, m2, n2],
            exponent: b,
            ..
        }: &GaussianPrimitive,
        diff: &Vector3<f64>,
        prod_center: &Vector3<f64>,
        point_charge: &PointCharge,
    ) -> f64 {
        let p = a + b;

        let diff_comp_nucleus = point_charge.position - prod_center;
        let dist_sq_comp_nucleus = diff_comp_nucleus.norm_squared();

        -point_charge.charge * std::f64::consts::TAU / p
            * (0..=l1 + l2)
                .flat_map(|t| {
                    let e1 = hermite_expansion([l1, l2, t], diff.x, a, b);
                    (0..=m1 + m2).flat_map(move |u| {
                        let e2 = hermite_expansion([m1, m2, u], diff.y, a, b);
                        (0..=n1 + n2).map(move |v| {
                            let e3 = hermite_expansion([n1, n2, v], diff.z, a, b);
                            e1 * e2
                                * e3
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

    /// four-electron integrals between four primitive gaussians
    fn _electron_repulsion(
        _primitives @ &[GaussianPrimitive {
            angular: [l1, m1, n1],
            exponent: a,
            ..
        }, GaussianPrimitive {
            angular: [l2, m2, n2],
            exponent: b,
            ..
        }, GaussianPrimitive {
            angular: [l3, m3, n3],
            exponent: c,
            ..
        }, GaussianPrimitive {
            angular: [l4, m4, n4],
            exponent: d,
            ..
        }]: &[GaussianPrimitive; 4],
        diff_ab: &Vector3<f64>,
        diff_cd: &Vector3<f64>,
        comp_diff: &Vector3<f64>,
        comp_dist_sq: f64,
    ) -> f64 {
        let p = a + b;
        let q = c + d;

        // TODO: make this more readable
        let alpha = p * q / (p + q);
        2.0 * std::f64::consts::PI.powi(5).sqrt()
            * (p * q * (p + q).sqrt()).recip()
            * (0..=l1 + l2)
                .flat_map(move |t1| {
                    let e1 = hermite_expansion([l1, l2, t1], diff_ab.x, a, b);
                    (0..=m1 + m2).flat_map(move |u1| {
                        let e2 = hermite_expansion([m1, m2, u1], diff_ab.y, a, b);
                        (0..=n1 + n2).flat_map(move |v1| {
                            let e3 = hermite_expansion([n1, n2, v1], diff_ab.z, a, b);
                            (0..=l3 + l4).flat_map(move |t2| {
                                let e4 = hermite_expansion([l3, l4, t2], diff_cd.x, c, d);
                                (0..=m3 + m4).flat_map(move |u2| {
                                    let e5 = hermite_expansion([m3, m4, u2], diff_cd.y, c, d);
                                    (0..=n3 + n4).map(move |v2| {
                                        let e6 = hermite_expansion([n3, n4, v2], diff_cd.z, c, d);
                                        e1 * e2
                                            * e3
                                            * e4
                                            * e5
                                            * e6
                                            * coulomb_auxiliary(
                                                t1 + t2,
                                                u1 + u2,
                                                v1 + v2,
                                                0,
                                                alpha,
                                                comp_diff,
                                                comp_dist_sq,
                                            )
                                            * if (t2 + u2 + v2) % 2 == 0 { 1.0 } else { -1.0 }
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
        let diff = b.position() - a.position();
        (0..a.primitives.len())
            .flat_map(move |i| {
                (0..b.primitives.len()).map(move |j| {
                    a.primitives[i].coefficient()
                        * b.primitives[j].coefficient()
                        * GaussianPrimitive::_overlap(&a.primitives[i], &b.primitives[j], &diff)
                })
            })
            .sum::<f64>()
    }

    fn kinetic_int(a: &Self, b: &Self) -> f64 {
        let diff = b.position() - a.position();
        (0..a.primitives.len())
            .flat_map(move |i| {
                (0..b.primitives.len()).map(move |j| {
                    a.primitives[i].coefficient()
                        * b.primitives[j].coefficient()
                        * GaussianPrimitive::_kinetic(&a.primitives[i], &b.primitives[j], &diff)
                })
            })
            .sum::<f64>()
    }

    fn nuclear_attraction_int(a: &Self, b: &Self, nuclei: &[PointCharge]) -> f64 {
        let m = a.primitives.len();
        let n = b.primitives.len();
        let diff = b.position() - a.position();
        nuclei
            .iter()
            .flat_map(move |point_charge| {
                (0..m).flat_map(move |i| {
                    (0..n).map(move |j| {
                        a.primitives[i].coefficient()
                            * b.primitives[j].coefficient()
                            * GaussianPrimitive::_nuclear_attraction(
                                &a.primitives[i],
                                &b.primitives[j],
                                &diff,
                                &GaussianPrimitive::product_center(
                                    &a.primitives[i],
                                    &a.position(),
                                    &b.primitives[j],
                                    &b.position(),
                                ),
                                point_charge,
                            )
                    })
                })
            })
            .sum::<f64>()
    }

    fn electron_repulsion_int(a: &Self, b: &Self, c: &Self, d: &Self) -> f64 {
        let a_pos = a.position();
        let b_pos = b.position();
        let c_pos = c.position();
        let d_pos = d.position();

        let diff_ab = b_pos - a_pos;
        let diff_cd = d_pos - c_pos;

        let itr = |pr: &Self| (0..pr.primitives.len());

        itr(a)
            .flat_map(move |i| {
                itr(b).flat_map(move |j| {
                    let comp_ab = GaussianPrimitive::product_center(
                        &a.primitives[i],
                        &a_pos,
                        &b.primitives[j],
                        &b_pos,
                    );

                    itr(c).flat_map(move |k| {
                        itr(d).map(move |l| {
                            let comp_cd = GaussianPrimitive::product_center(
                                &c.primitives[k],
                                &c_pos,
                                &d.primitives[l],
                                &d_pos,
                            );

                            let comp_diff = comp_cd - comp_ab;
                            let comp_dist_sq = comp_diff.norm_squared();

                            a.primitives[i].coefficient()
                                * b.primitives[j].coefficient()
                                * c.primitives[k].coefficient()
                                * d.primitives[l].coefficient()
                                * GaussianPrimitive::_electron_repulsion(
                                    &[
                                        a.primitives[i],
                                        b.primitives[j],
                                        c.primitives[k],
                                        d.primitives[l],
                                    ],
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
