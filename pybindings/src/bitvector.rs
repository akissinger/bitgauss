use pyo3::exceptions::PyValueError;
use pyo3::{prelude::*, IntoPyObjectExt};

use bitgauss::bitvector::BitVector;
use bitgauss::BitMatrix;
use rand::{rngs::SmallRng, SeedableRng};

use crate::bitmatrix::PyBitMatrix;

#[pyclass(name = "BitVector")]
#[derive(Clone)]
pub struct PyBitVector {
    pub(crate) inner: BitVector,
}

#[pymethods]
impl PyBitVector {
    /// Creates a new BitVector of specified length initialized to zero
    #[new]
    pub fn new(length: usize) -> Self {
        PyBitVector {
            inner: BitVector::zeros(length),
        }
    }

    /// Gets the bit at position i
    pub fn bit(&self, i: usize) -> PyResult<bool> {
        if i >= self.inner.len() {
            return Err(PyValueError::new_err("Index out of bounds"));
        }
        Ok(self.inner.bit(i))
    }

    /// Sets the bit at position i to b
    pub fn set_bit(&mut self, i: usize, b: bool) -> PyResult<()> {
        if i >= self.inner.len() {
            return Err(PyValueError::new_err("Index out of bounds"));
        }
        self.inner.set_bit(i, b);
        Ok(())
    }

    /// Builds a BitVector from a Python function that determines the value of each bit
    #[staticmethod]
    pub fn build(length: usize, func: PyObject) -> PyResult<Self> {
        Python::with_gil(|py| {
            let vector = BitVector::build(length, |i| {
                let result = func.call1(py, (i,));
                match result {
                    Ok(val) => val.is_truthy(py).unwrap_or(false),
                    Err(_) => false,
                }
            });
            Ok(PyBitVector { inner: vector })
        })
    }

    /// Creates a new BitVector of specified length initialized to zero
    #[staticmethod]
    pub fn zeros(length: usize) -> Self {
        PyBitVector {
            inner: BitVector::zeros(length),
        }
    }

    /// Checks if the vector consists of all zero bits
    pub fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    /// Creates a new random BitVector of specified length
    #[staticmethod]
    #[pyo3(signature = (length, seed=None))]
    pub fn random(length: usize, seed: Option<u64>) -> Self {
        let mut rng = if let Some(s) = seed {
            SmallRng::seed_from_u64(s)
        } else {
            SmallRng::from_os_rng()
        };

        PyBitVector {
            inner: BitVector::random(&mut rng, length),
        }
    }

    /// Returns the length of the vector
    #[getter]
    pub fn length(&self) -> usize {
        self.inner.len()
    }

    /// Returns true if the vector has length 0
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Returns the number of 1s in the vector (Hamming weight)
    pub fn weight(&self) -> usize {
        self.inner.weight()
    }

    /// XORs another BitVector into this one
    pub fn xor_with(&mut self, other: &PyBitVector) -> PyResult<()> {
        if self.inner.len() != other.inner.len() {
            return Err(PyValueError::new_err(
                "BitVectors must have the same length for XOR",
            ));
        }
        self.inner.xor_with(&other.inner);
        Ok(())
    }

    /// Returns a copy of the vector
    pub fn copy(&self) -> Self {
        PyBitVector {
            inner: self.inner.clone(),
        }
    }

    /// String representation of the vector
    pub fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Python representation of the vector
    pub fn __repr__(&self) -> String {
        format!("BitVector(length={})", self.inner.len())
    }

    /// Support for indexing with [i]
    pub fn __getitem__(&self, key: PyObject) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            if let Ok(i) = key.extract::<usize>(py) {
                // Vector[i] - get single bit
                if i >= self.inner.len() {
                    return Err(PyValueError::new_err("Index out of bounds"));
                }
                self.inner.bit(i).into_py_any(py)
            } else {
                Err(PyValueError::new_err("Invalid index type"))
            }
        })
    }

    /// Support for item assignment with [i] = value
    pub fn __setitem__(&mut self, key: PyObject, value: PyObject) -> PyResult<()> {
        Python::with_gil(|py| {
            if let Ok(i) = key.extract::<usize>(py) {
                if i >= self.inner.len() {
                    return Err(PyValueError::new_err("Index out of bounds"));
                }
                let bit_value = value.is_truthy(py)?;
                self.inner.set_bit(i, bit_value);
                Ok(())
            } else {
                Err(PyValueError::new_err("Invalid index type for assignment"))
            }
        })
    }

    /// XOR operation using the ^ operator
    pub fn __xor__(&self, other: &PyBitVector) -> PyResult<Self> {
        if self.inner.len() != other.inner.len() {
            return Err(PyValueError::new_err(
                "BitVectors must have the same length for XOR",
            ));
        }
        let result = &self.inner ^ &other.inner;
        Ok(PyBitVector { inner: result })
    }

    /// Right-hand XOR operation
    pub fn __rxor__(&self, other: &PyBitVector) -> PyResult<Self> {
        other.__xor__(self)
    }

    /// In-place XOR operation using ^=
    pub fn __ixor__(&mut self, other: &PyBitVector) -> PyResult<()> {
        self.xor_with(other)
    }

    /// Vector equality comparison
    pub fn __eq__(&self, other: &PyBitVector) -> bool {
        self.inner == other.inner
    }

    /// Vector inequality comparison
    pub fn __ne__(&self, other: &PyBitVector) -> bool {
        !self.__eq__(other)
    }

    /// Returns the length of the vector (for len() function)
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Convert vector to a list of bools
    pub fn to_list(&self) -> Vec<bool> {
        (0..self.inner.len()).map(|i| self.inner.bit(i)).collect()
    }

    /// Create vector from a list of bools
    #[staticmethod]
    pub fn from_list(data: Vec<bool>) -> Self {
        PyBitVector {
            inner: BitVector::from_bool_vec(&data),
        }
    }

    /// Convert vector to a list of integers (0 or 1)
    pub fn to_int_list(&self) -> Vec<usize> {
        (0..self.inner.len())
            .map(|i| if self.inner.bit(i) { 1 } else { 0 })
            .collect()
    }

    /// Create vector from a list of integers (0 or 1)
    #[staticmethod]
    pub fn from_int_list(data: Vec<usize>) -> Self {
        PyBitVector {
            inner: BitVector::from_int_vec(&data),
        }
    }

    /// Convert to BitMatrix (single row)
    pub fn to_matrix(&self) -> PyBitMatrix {
        PyBitMatrix::from(BitMatrix::from(self.inner.clone()))
    }

    /// Create BitVector from a single-row BitMatrix
    #[staticmethod]
    pub fn from_matrix(matrix: &PyBitMatrix) -> PyResult<Self> {
        let bit_matrix = &matrix.inner;
        match BitVector::try_from(bit_matrix.clone()) {
            Ok(vector) => Ok(PyBitVector { inner: vector }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Cannot convert BitMatrix to BitVector: {}",
                e
            ))),
        }
    }
}

// Helper implementations

impl From<BitVector> for PyBitVector {
    fn from(inner: BitVector) -> Self {
        PyBitVector { inner }
    }
}

impl From<PyBitVector> for BitVector {
    fn from(py_vector: PyBitVector) -> Self {
        py_vector.inner
    }
}
