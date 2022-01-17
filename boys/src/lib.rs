use crate::lut::LUT;

pub mod lut;

fn double_factorial(n: u32) -> u32 {
    match n {
        n @ 1..=u32::MAX => (2 - n % 2..=n).step_by(2).product(),
        _ => 1,
    }
}

// see: https://sci-hub.mksa.top/https://doi.org/10.1002/jcc.23935
fn _boys_impl(n: i32, t: f64, mut curr: f64, exp_mt: f64) -> f64 {
    // downwards recursion
    let mut curr_n = MAX_TABULATED_N;
    while curr_n > n {
        curr = (2.0 * t * curr + exp_mt) / (2 * curr_n - 1) as f64;
        curr_n -= 1;
    }
    curr
}

const DT: f64 = 0.15;
const MAX_TABULATED_T: f64 = 30.0;
const MAX_TABULATED_N: i32 = 10;
pub(crate) const N: usize = (MAX_TABULATED_T as f64 / DT) as usize;

pub fn boys_quadrature(t: f64, n: i32) -> f64 {
    let integrand = |x: f64| -> f64 { x.powi(2 * n) * f64::exp(-t * x.powi(2)) };

    const ITERS: i32 = 15;
    const ITERS_F: f64 = ITERS as f64;

    ITERS_F.recip()
        * (integrand(0.0) * 0.5
            + (1..ITERS)
                .map(|k| integrand(k as f64 * ITERS_F.recip()))
                .sum::<f64>()
            + integrand(1.0) * 0.5)
}

pub fn boys(t: f64, n: i32) -> f64 {
    assert!(
        n <= MAX_TABULATED_N,
        "Boys-Function order {} is higher than max tabulated order {}",
        n,
        MAX_TABULATED_N
    );

    if t.abs() < 1e-16 {
        f64::recip((2 * n + 1) as f64)
    } else if t.abs() >= MAX_TABULATED_T {
        double_factorial((2 * n as u32).saturating_sub(1)) as f64
            * f64::powi(2.0, n + 1).recip()
            * f64::sqrt(std::f64::consts::PI * t.powi(-2 * n - 1))
    } else {
        let index_f = N as f64 * t / MAX_TABULATED_T;
        let lower = index_f.floor();
        let upper = index_f.ceil();
        let int_t = index_f - lower;

        let lower = lower as usize;
        let upper = upper as usize;
        let curr_val = int_t * LUT[upper.min(N - 1)] + (1.0 - int_t) * LUT[lower];
        _boys_impl(n, t, curr_val, f64::exp(-t))
    }
}

#[test]
pub fn test_boys() {
    use quadrature::integrate;
    for n in 0..=MAX_TABULATED_N {
        let mut max_diff = 0.0f64;
        let mut avg_diff = 0.0f64;
        for t in (0..35_000).map(|k| k as f64 / 1000.0) {
            let result = integrate(
                |x| x.powi(2 * n) * f64::exp(-t * x.powi(2)),
                0.0,
                1.0,
                1e-32,
            );
            let real = result.integral;
            let acq = boys(t, n);
            let diff = (real - acq).abs();
            max_diff = max_diff.max(diff);
            avg_diff += diff;
        }

        println!(
            "Degree {}: mean diff {:0.5e}, max diff {:0.5e}",
            n,
            avg_diff / 35_000.0,
            max_diff
        );
    }
}

#[test]
pub fn gen_lut() {
    use quadrature::integrate;
    const EPSILON: f64 = 1e-64;

    print!("pub(crate) const LUT: &[f64; crate::N] = &[");
    for t in (0..N).map(|n| n as f64 * DT) {
        let result = integrate(
            |x| x.powi(2 * MAX_TABULATED_N) * f64::exp(-t * x.powi(2)),
            0.0,
            1.0,
            EPSILON,
        );
        print!("{:e},", result.integral);
    }
    print!("];");
}
