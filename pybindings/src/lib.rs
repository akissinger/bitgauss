// There seems to be some issues with the pyo3 bindings generation on methods returning
// a `PyResult<T>`.
#![allow(
    clippy::useless_conversion,
    clippy::uninlined_format_args,
    clippy::must_use_candidate,
    clippy::return_self_not_must_use,
    clippy::needless_pass_by_value,
    clippy::missing_errors_doc,
    clippy::bool_to_int_with_if
)]

pub mod bitmatrix;
pub mod bitvector;

use crate::bitmatrix::PyBitMatrix;
use crate::bitvector::PyBitVector;
use pyo3::prelude::*;

#[pymodule]
fn bitgauss(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBitMatrix>()?;
    m.add_class::<PyBitVector>()?;
    Ok(())
}
