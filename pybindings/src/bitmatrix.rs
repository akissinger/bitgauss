use pyo3::prelude::*;

#[pyclass]
pub struct BitMatrix(pub(crate) bitgauss::BitMatrix);
