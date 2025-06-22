// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(clippy::useless_conversion)]

pub mod bitmatrix;

use crate::bitmatrix::PyBitMatrix;
use pyo3::prelude::*;

#[pymodule]
fn bitgauss(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBitMatrix>()?;
    Ok(())
}
