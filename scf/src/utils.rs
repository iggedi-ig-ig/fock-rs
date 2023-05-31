use std::{collections::VecDeque, mem::MaybeUninit};

use faer_evd::SymmetricEvdParams;
use itertools::iproduct;
use nalgebra::{DMatrix, DVector, SymmetricEigen};

static mut BUFFER: [MaybeUninit<u8>; 2 << 18] = [MaybeUninit::uninit(); 2 << 18];

/// Constructs a Hermitian matrix of size `n` by evaluating the function `func` for all pairs of indices
/// `(i, j)` where `0 <= i < n` and `0 <= j < n`.
///
/// The function `func` is only called for indices where `i <= j`. The resulting upper-triangular matrix is
/// then mirrored across the diagonal to produce a full Hermitian matrix.
///
/// # Arguments
///
/// * `n` - The size of the matrix to construct.
/// * `func` - A closure that takes two `usize` arguments `i` and `j`, and returns a `f64` value that represents
///            the entry at row `i` and column `j` in the resulting Hermitian matrix.
///
/// # Returns
///
/// A `DMatrix<f64>` object representing a Hermitian matrix of size `n` that was constructed by evaluating
/// the function `func` for all pairs of indices `(i, j)` where `0 <= i < n` and `0 <= j < n`.
pub fn hermitian(n: usize, mut func: impl FnMut(usize, usize) -> f64) -> DMatrix<f64> {
    let m = DMatrix::from_fn(n, n, |i, j| if i <= j { func(i, j) } else { 0.0 });
    DMatrix::from_fn(n, n, |i, j| if i <= j { m[(i, j)] } else { m[(j, i)] })
}

/// Computes the eigenvectors and eigenvalues of a given real symmetric matrix using the `SymmetricEigen` type of nalgebra.
///
/// This function is a wrapper around `SymmetricEigen`, and requires the input matrix to be real symmetric for accurate results.
///
/// # Arguments
///
/// * `matrix` - The `DMatrix<f64>` representing the real symmetric matrix for which to compute the eigenvectors and eigenvalues.
///
/// # Returns
///
/// A tuple `(eigenvectors, eigenvalues)` where `eigenvectors` is a `DMatrix<f64>` containing the eigenvectors of the input matrix,
/// and `eigenvalues` is a `DVector<f64>` containing the eigenvalues of the input matrix.
pub fn eigs(matrix: DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
    let eigs = SymmetricEigen::new(matrix);
    (eigs.eigenvectors, eigs.eigenvalues)
}

/// Computes the eigenvectors and eigenvalues of a given real symmetric matrix using the `SymmetricEigen` type of nalgebra,
/// and sorts the eigenvector - eigenvalue pairs by the eigenvalues in an ascending order.
///
/// # Arguments
///
/// * `matrix` - The `DMatrix<f64>` representing the real symmetric matrix for which to compute the eigenvectors and eigenvalues.
///
/// # Returns
///
/// A tuple `(eigenvectors, eigenvalues)` where `eigenvectors` is a `DMatrix<f64>` containing the eigenvectors of the input matrix,
/// sorted by the eigenvalues in an ascending order, and `eigenvalues` is a `DVector<f64>` containing the eigenvalues of the input matrix,
/// also sorted in ascending order.
pub fn sorted_eigs(matrix: DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
    let nrows = matrix.nrows();
    let ncols = matrix.ncols();

    let stack = dyn_stack::DynStack::new(unsafe { &mut BUFFER });
    let matrix = faer_core::Mat::with_dims(nrows, ncols, |i, j| matrix[(i, j)]);
    let mut s = faer_core::Mat::zeros(nrows, 1);
    let mut u = faer_core::Mat::zeros(nrows, ncols);

    faer_evd::compute_hermitian_evd(
        matrix.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        1e-8,
        1e-8,
        faer_core::Parallelism::Rayon(4),
        stack,
        SymmetricEvdParams::default(),
    );

    let eigenvectors = DMatrix::from_columns(
        &(0..ncols)
            .map(|j| DVector::from_iterator(nrows, (0..nrows).map(|i| u.read(i, j))))
            .collect::<Vec<_>>(),
    );
    let eigenvalues = DVector::from_iterator(ncols, (0..ncols).map(|i| s.read(i, 0)));

    let mut val_vec_pairs = eigenvalues
        .into_iter()
        .zip(eigenvectors.column_iter())
        .collect::<Vec<_>>();

    val_vec_pairs.sort_unstable_by(|(a, _), (b, _)| a.total_cmp(b));

    let (values, vectors): (Vec<_>, Vec<_>) = val_vec_pairs.into_iter().unzip();

    (
        DMatrix::from_columns(&vectors),
        DVector::from_column_slice(&values),
    )
}

/// Performs DIIS (direct inversion of iterative subspace), given a list of error vectors and fock matricies and returns a new fock matrix.
pub fn diis(
    error_vectors: &VecDeque<DMatrix<f64>>,
    fock_matricies: &VecDeque<DMatrix<f64>>,
) -> Option<DMatrix<f64>> {
    assert_eq!(error_vectors.len(), fock_matricies.len());
    let n = error_vectors.len();

    let mut matrix = DMatrix::zeros(n + 1, n + 1);
    // upper block
    for (i, j) in iproduct!(0..n, 0..n) {
        matrix[(i, j)] = error_vectors[i].dot(&error_vectors[j]);
    }

    // last row
    for i in 0..n {
        matrix[(n, i)] = -1.0;
    }

    // last col
    for i in 0..n {
        matrix[(i, n)] = -1.0;
    }

    // last entry
    matrix[(n, n)] = 0.0;

    let mut b = DVector::zeros(n + 1);
    b[(n, 0)] = -1.0;

    matrix.try_inverse().map(|inv| inv * b).map(|c| {
        c.iter()
            .enumerate()
            .take(n)
            .map(|(i, &x)| x * &fock_matricies[i])
            .sum()
    })
}
