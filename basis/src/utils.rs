use nalgebra::Vector3;

#[allow(clippy::many_single_char_names)]
// TODO: make this iterative
// see https://joshuagoings.com/2017/04/28/integrals/
pub fn hermite_expansion([i, j, t]: [i32; 3], dist: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    let q = a * b / p;

    if t < 0 || t > i + j {
        0.0
    } else if i == j && j == t && t == 0 {
        f64::exp(-q * dist.powi(2))
    } else if j == 0 {
        (2.0 * p).recip() * hermite_expansion([i - 1, j, t - 1], dist, a, b)
            - (q * dist / a) * hermite_expansion([i - 1, j, t], dist, a, b)
            + (t + 1) as f64 * hermite_expansion([i - 1, j, t + 1], dist, a, b)
    } else {
        (2.0 * p).recip() * hermite_expansion([i, j - 1, t - 1], dist, a, b)
            + (q * dist / b) * hermite_expansion([i, j - 1, t], dist, a, b)
            + (t + 1) as f64 * hermite_expansion([i, j - 1, t + 1], dist, a, b)
    }
}

#[allow(clippy::many_single_char_names)]
// TODO: make this iterative
// see https://joshuagoings.com/2017/04/28/integrals/
pub fn coulomb_auxiliary(
    t: i32,
    u: i32,
    v: i32,
    n: i32,
    p: f64,
    diff: &Vector3<f64>,
    dist_sq: f64,
) -> f64 {
    if t == u && u == v && v == 0 {
        (-2.0 * p).powi(n) * boys::boys(p * dist_sq, n)
    } else if t == u && u == 0 {
        diff.z * coulomb_auxiliary(t, u, v - 1, n + 1, p, diff, dist_sq)
            + if v > 1 {
                (v - 1) as f64 * coulomb_auxiliary(t, u, v - 2, n + 1, p, diff, dist_sq)
            } else {
                0.0
            }
    } else if t == 0 {
        diff.y * coulomb_auxiliary(t, u - 1, v, n + 1, p, diff, dist_sq)
            + if u > 1 {
                (u - 1) as f64 * coulomb_auxiliary(t, u - 2, v, n + 1, p, diff, dist_sq)
            } else {
                0.0
            }
    } else {
        diff.x * coulomb_auxiliary(t - 1, u, v, n + 1, p, diff, dist_sq)
            + if t > 1 {
                coulomb_auxiliary(t - 2, u, v, n + 1, p, diff, dist_sq)
            } else {
                0.0
            }
    }
}
