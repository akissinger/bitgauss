use pyo3::exceptions::PyValueError;
use pyo3::{prelude::*, IntoPyObjectExt};

// Assuming these are the original imports/dependencies from your Rust code
// You'll need to adjust these based on your actual crate structure
use bitgauss::BitMatrix;
use rand::{rngs::SmallRng, SeedableRng};

#[pyclass(name = "BitMatrix")]
pub struct PyBitMatrix {
    inner: BitMatrix,
}

#[pymethods]
impl PyBitMatrix {
    /// Creates a new BitMatrix of size rows x cols with all bits set to 0
    #[new]
    pub fn new(rows: usize, cols: usize) -> Self {
        PyBitMatrix {
            inner: BitMatrix::zeros(rows, cols),
        }
    }

    /// Gets the bit at position (i, j)
    pub fn bit(&self, i: usize, j: usize) -> PyResult<bool> {
        if i >= self.inner.rows() || j >= self.inner.cols() {
            return Err(PyValueError::new_err("Index out of bounds"));
        }
        Ok(self.inner.bit(i, j))
    }

    /// Sets the bit at position (i, j) to b
    pub fn set_bit(&mut self, i: usize, j: usize, b: bool) -> PyResult<()> {
        if i >= self.inner.rows() || j >= self.inner.cols() {
            return Err(PyValueError::new_err("Index out of bounds"));
        }
        self.inner.set_bit(i, j, b);
        Ok(())
    }

    /// Builds a BitMatrix from a Python function that determines the value of each bit
    #[staticmethod]
    pub fn build(rows: usize, cols: usize, func: PyObject) -> PyResult<Self> {
        Python::with_gil(|py| {
            let matrix = BitMatrix::build(rows, cols, |i, j| {
                let result = func.call1(py, (i, j));
                match result {
                    Ok(val) => val.is_truthy(py).unwrap_or(false),
                    Err(_) => false,
                }
            });
            Ok(PyBitMatrix { inner: matrix })
        })
    }

    /// Creates a new BitMatrix of size rows x cols with all bits set to 0
    #[staticmethod]
    pub fn zeros(rows: usize, cols: usize) -> Self {
        PyBitMatrix {
            inner: BitMatrix::zeros(rows, cols),
        }
    }

    /// Checks if the matrix consists of all zero bits
    pub fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    /// Creates a new identity BitMatrix of size size x size
    #[staticmethod]
    pub fn identity(size: usize) -> Self {
        PyBitMatrix {
            inner: BitMatrix::identity(size),
        }
    }

    /// Creates a new random BitMatrix of size rows x cols
    #[staticmethod]
    #[pyo3(signature = (rows, cols, seed=None))]
    pub fn random(rows: usize, cols: usize, seed: Option<u64>) -> Self {
        let mut rng = if let Some(s) = seed {
            SmallRng::seed_from_u64(s)
        } else {
            SmallRng::from_os_rng()
        };

        PyBitMatrix {
            inner: BitMatrix::random(&mut rng, rows, cols),
        }
    }

    /// Creates a new random invertible BitMatrix of size size x size
    #[staticmethod]
    #[pyo3(signature = (size, seed=None))]
    pub fn random_invertible(size: usize, seed: Option<u64>) -> Self {
        let mut rng = if let Some(s) = seed {
            SmallRng::seed_from_u64(s)
        } else {
            SmallRng::from_os_rng()
        };

        PyBitMatrix {
            inner: BitMatrix::random_invertible(&mut rng, size),
        }
    }

    /// Returns the number of logical rows in the matrix
    #[getter]
    pub fn rows(&self) -> usize {
        self.inner.rows()
    }

    /// Returns the number of logical columns in the matrix
    #[getter]
    pub fn cols(&self) -> usize {
        self.inner.cols()
    }

    /// Returns a transposed copy of the matrix
    pub fn transposed(&self) -> Self {
        PyBitMatrix {
            inner: self.inner.transposed(),
        }
    }

    /// Transposes the matrix in place
    pub fn transpose_inplace(&mut self) {
        self.inner.transpose_inplace();
    }

    /// Performs gaussian elimination
    #[pyo3(signature = (full=false))]
    pub fn gauss(&mut self, full: bool) {
        self.inner.gauss(full);
    }

    /// Computes the rank of the matrix using gaussian elimination
    pub fn rank(&self) -> usize {
        self.inner.rank()
    }

    /// Computes the inverse of an invertible matrix
    pub fn inverse(&self) -> PyResult<Self> {
        match self.inner.try_inverse() {
            Ok(inv) => Ok(PyBitMatrix { inner: inv }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Matrix inversion failed: {}",
                e
            ))),
        }
    }

    /// Vertically stacks this matrix with another one and returns the result
    pub fn vstack(&self, other: &PyBitMatrix) -> PyResult<Self> {
        match self.inner.try_vstack(&other.inner) {
            Ok(result) => Ok(PyBitMatrix { inner: result }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Vertical stack failed: {}",
                e
            ))),
        }
    }

    /// Horizontally stacks this matrix with another one and returns the result
    pub fn hstack(&self, other: &PyBitMatrix) -> PyResult<Self> {
        match self.inner.try_hstack(&other.inner) {
            Ok(result) => Ok(PyBitMatrix { inner: result }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Horizontal stack failed: {}",
                e
            ))),
        }
    }

    /// Computes a basis for the nullspace of the matrix
    pub fn nullspace(&self) -> Vec<PyBitMatrix> {
        self.inner
            .nullspace()
            .into_iter()
            .map(|m| PyBitMatrix { inner: m })
            .collect()
    }

    /// Returns a copy of the matrix
    pub fn copy(&self) -> Self {
        PyBitMatrix {
            inner: self.inner.clone(),
        }
    }

    /// String representation of the matrix
    pub fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Python representation of the matrix
    pub fn __repr__(&self) -> String {
        format!("BitMatrix({}x{})", self.inner.rows(), self.inner.cols())
    }

    /// Support for indexing with [i, j]
    pub fn __getitem__(&self, key: PyObject) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            if let Ok((i, j)) = key.extract::<(usize, usize)>(py) {
                // Matrix[i, j] - get single bit
                if i >= self.inner.rows() || j >= self.inner.cols() {
                    return Err(PyValueError::new_err("Index out of bounds"));
                }
                self.inner.bit(i, j).into_py_any(py)
            } else {
                Err(PyValueError::new_err("Invalid index type"))
            }
        })
    }

    /// Support for item assignment with [i, j] = value
    pub fn __setitem__(&mut self, key: PyObject, value: PyObject) -> PyResult<()> {
        Python::with_gil(|py| {
            if let Ok((i, j)) = key.extract::<(usize, usize)>(py) {
                if i >= self.inner.rows() || j >= self.inner.cols() {
                    return Err(PyValueError::new_err("Index out of bounds"));
                }
                let bit_value = value.is_truthy(py)?;
                self.inner.set_bit(i, j, bit_value);
                Ok(())
            } else {
                Err(PyValueError::new_err("Invalid index type for assignment"))
            }
        })
    }

    /// Matrix multiplication using the @ operator (Python 3.5+)
    pub fn __matmul__(&self, other: &PyBitMatrix) -> PyResult<Self> {
        match self.inner.try_mul(&other.inner) {
            Ok(result) => Ok(PyBitMatrix { inner: result }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Matrix multiplication failed: {}",
                e
            ))),
        }
    }

    /// Matrix multiplication using the * operator (for compatibility)
    pub fn __mul__(&self, other: &PyBitMatrix) -> PyResult<Self> {
        self.__matmul__(other)
    }

    /// Right-hand matrix multiplication
    pub fn __rmul__(&self, other: &PyBitMatrix) -> PyResult<Self> {
        other.__matmul__(self)
    }

    /// In-place matrix multiplication using @=
    pub fn __imatmul__(&mut self, other: &PyBitMatrix) -> PyResult<()> {
        let result = self.__matmul__(other)?;
        self.inner = result.inner;
        Ok(())
    }

    /// Matrix multiplication method (alternative to operators)
    pub fn matmul(&self, other: &PyBitMatrix) -> PyResult<Self> {
        self.__matmul__(other)
    }

    /// Matrix power (repeated multiplication)
    pub fn __pow__(&self, exponent: usize, _modulus: Option<&PyBitMatrix>) -> PyResult<Self> {
        if self.inner.rows() != self.inner.cols() {
            return Err(PyValueError::new_err(
                "Matrix power is only defined for square matrices",
            ));
        }

        if exponent == 0 {
            return Ok(PyBitMatrix {
                inner: BitMatrix::identity(self.inner.rows()),
            });
        }

        let mut result = self.copy();
        for _ in 1..exponent {
            result.inner = &result.inner * &self.inner;
        }

        Ok(result)
    }

    /// Convert matrix to a list of lists of bools
    pub fn to_list(&self) -> Vec<Vec<bool>> {
        (0..self.inner.rows())
            .map(|i| {
                (0..self.inner.cols())
                    .map(|j| self.inner.bit(i, j))
                    .collect()
            })
            .collect()
    }

    /// Create matrix from a list of lists of bools
    #[staticmethod]
    pub fn from_list(data: Vec<Vec<bool>>) -> PyResult<Self> {
        if data.is_empty() {
            return Ok(PyBitMatrix::zeros(0, 0));
        }

        let rows = data.len();
        let cols = data[0].len();

        // Check that all rows have the same length
        for (i, row) in data.iter().enumerate() {
            if row.len() != cols {
                return Err(PyValueError::new_err(format!(
                    "Row {} has length {}, expected {}",
                    i,
                    row.len(),
                    cols
                )));
            }
        }

        let matrix = BitMatrix::build(rows, cols, |i, j| data[i][j]);
        Ok(PyBitMatrix { inner: matrix })
    }

    /// Convert matrix to a list of lists of integers (0 or 1)
    pub fn to_int_list(&self) -> Vec<Vec<usize>> {
        (0..self.inner.rows())
            .map(|i| {
                (0..self.inner.cols())
                    .map(|j| if self.inner.bit(i, j) { 1 } else { 0 })
                    .collect()
            })
            .collect()
    }

    /// Create matrix from a list of lists of integers (0 or 1)
    #[staticmethod]
    pub fn from_int_list(data: Vec<Vec<usize>>) -> PyResult<Self> {
        if data.is_empty() {
            return Ok(PyBitMatrix::zeros(0, 0));
        }

        let rows = data.len();
        let cols = data[0].len();

        // Check that all rows have the same length
        for (i, row) in data.iter().enumerate() {
            if row.len() != cols {
                return Err(PyValueError::new_err(format!(
                    "Row {} has length {}, expected {}",
                    i,
                    row.len(),
                    cols
                )));
            }
        }

        let matrix = BitMatrix::build(rows, cols, |i, j| data[i][j] != 0);
        Ok(PyBitMatrix { inner: matrix })
    }

    /// Matrix equality comparison
    pub fn __eq__(&self, other: &PyBitMatrix) -> bool {
        self.inner == other.inner
    }

    /// Matrix inequality comparison
    pub fn __ne__(&self, other: &PyBitMatrix) -> bool {
        !self.__eq__(other)
    }
}

// Additional helper implementations that might be needed

impl From<BitMatrix> for PyBitMatrix {
    fn from(inner: BitMatrix) -> Self {
        PyBitMatrix { inner }
    }
}

impl From<PyBitMatrix> for BitMatrix {
    fn from(py_matrix: PyBitMatrix) -> Self {
        py_matrix.inner
    }
}
