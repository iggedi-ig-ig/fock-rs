use nalgebra::{DMatrix, DVector, SymmetricEigen};
use std::cmp::Ordering;

/// constructs a hermitian matrix. This only calculates the upper triangular parts and then mirrors it
pub fn hermitian(n: usize, mut func: impl FnMut(usize, usize) -> f64) -> DMatrix<f64> {
    let m = DMatrix::from_fn(n, n, |i, j| if i <= j { func(i, j) } else { 0.0 });
    DMatrix::from_fn(n, n, |i, j| if i <= j { m[(i, j)] } else { m[(j, i)] })
}

/// returns eigenvector - eigenvalue pairs of the given matrix
pub fn eigs(matrix: DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
    let eigs = SymmetricEigen::new(matrix);
    (eigs.eigenvectors, eigs.eigenvalues)
}

/// returns eigenvector - eigenvalue pairs, sorted by the eigenvalues in an ascending order
pub fn sorted_eigs(matrix: DMatrix<f64>) -> (DMatrix<f64>, DVector<f64>) {
    let (eigenvectors, eigenvalues) = eigs(matrix);

    let mut val_vec_pairs = eigenvalues
        .into_iter()
        .zip(eigenvectors.column_iter())
        .collect::<Vec<_>>();

    val_vec_pairs.sort_unstable_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(Ordering::Less));

    let (values, vectors): (Vec<_>, Vec<_>) = val_vec_pairs.into_iter().unzip();

    (
        DMatrix::from_columns(&vectors),
        DVector::from_column_slice(&values),
    )
}
