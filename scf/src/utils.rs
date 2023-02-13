use nalgebra::{DMatrix, DVector, SymmetricEigen};

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
    let (eigenvectors, eigenvalues) = eigs(matrix);

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
