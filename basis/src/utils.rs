use nalgebra::Vector3;

/// Calculates the Hermite expansion of two arbitrary Gaussian functions centered at a and b with exponents i and j,
/// respectively, for a distance of `dist`. The expansion is calculated at order `t`. This function is not optimized for
/// speed and should be improved to use iterative methods.
///
/// # Arguments
///
/// * `[i, j, t]` - An array of exponents and expansion order
/// * `dist` - The distance between the two Gaussian functions
/// * `a` - The exponent of the first Gaussian function
/// * `b` - The exponent of the second Gaussian function
///
/// # Returns
///
/// The value of the Hermite expansion at order `t`. If `t` is less than zero or greater
/// than `i + j`, this function will return `0.0`.
///
/// # Examples
///
/// ```
/// use basis::utils::hermite_expansion;
///
/// assert_eq!(hermite_expansion([1, 2, 3], 2.5, 3.0, 5.0), 1.9872458203327215e-9);
/// ```
///
/// # References
///
/// [1] Goings, J. Integrals. https://joshuagoings.com/2017/04/28/integrals/
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

/// coulomb_auxiliary calculates the auxiliary integral for two-electron integrals in Gaussian type
/// orbital basis functions (GTOs).
///
/// # Arguments
///
/// * `t` - an `i32` representing the exponent of the x-component of the first GTO.
/// * `u` - an `i32` representing the exponent of the y-component of the first GTO.
/// * `v` - an `i32` representing the exponent of the z-component of the first GTO.
/// * `n` - an `i32` representing the order of the auxiliary integral.
/// * `p` - a `f64` representing the sum of the two exponents of the contracted GTOs.
/// * `diff` - a `&Vector3<f64>` representing the difference vector between the center of the
///            contracted GTOs.
/// * `dist_sq` - a `f64` representing the squared distance between the centers of the contracted
///               GTOs.
///
/// # Returns
///
/// The value of the auxiliary integral.
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
