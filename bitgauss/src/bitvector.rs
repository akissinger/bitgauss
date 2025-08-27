use crate::bitmatrix::BitMatrix;
use std::fmt;
use std::ops::{BitXor, BitXorAssign, Index, Mul};

/// A wrapper around a one-row `BitMatrix`
///
/// Despite storing bits in row-major order, this struct behaves like a column vector for the
/// purposes of matrix multiplication. To multiply as a row vector instead, use
/// `BitVector::as_row_vector()`.
#[derive(Clone, Debug)]
pub struct BitVector(BitMatrix);

impl BitVector {
    /// Gets the bit at position `i`
    #[inline]
    pub fn bit(&self, i: usize) -> bool {
        self.0.bit(0, i)
    }

    /// Sets the bit at position `i` to `b`
    #[inline]
    pub fn set_bit(&mut self, i: usize, b: bool) {
        self.0.set_bit(0, i, b);
    }

    /// Builds a `BitVector` from a function `f` that determines the value of each bit
    ///
    /// # Arguments
    /// * `length` - the number of columns in the matrix
    /// * `f` - a function that takes the row and column indices and returns a boolean value for each bit
    pub fn build(length: usize, mut f: impl FnMut(usize) -> bool) -> Self {
        Self(BitMatrix::build(1, length, |_, j| f(j)))
    }

    /// Creates a new `BitVector` from a vector of bool vectors
    pub fn from_bool_vec(data: &Vec<bool>) -> Self {
        Self::build(data.len(), |i| data[i])
    }

    /// Creates a new `BitVector` from a vector of integer vectors
    pub fn from_int_vec(data: &Vec<usize>) -> Self {
        Self::build(data.len(), |i| data[i] != 0)
    }

    /// Creates a new `BitVector` of size `length` initialized to zero
    pub fn zeros(length: usize) -> Self {
        Self(BitMatrix::zeros(1, length))
    }

    /// Checks if the vector consists of all zero bits
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Returns the length of the vector
    #[inline]
    pub fn len(&self) -> usize {
        self.0.cols()
    }

    /// Returns true if the vector has length 0
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Creates a new random `BitVector` of specified length
    #[inline]
    pub fn random(rng: &mut impl rand::Rng, length: usize) -> Self {
        Self(BitMatrix::random(rng, 1, length))
    }

    /// Returns the number of 1s in the vector (Hamming weight)
    #[inline]
    pub fn weight(&self) -> usize {
        self.0.row_weight(0)
    }

    /// XORs another `BitVector` into this one
    #[inline]
    pub fn xor_with(&mut self, other: &BitVector) {
        if self.len() != other.len() {
            panic!("BitVectors must have the same length for XOR");
        }
        self.0.add_bits_to_row(other.0.row(0), 0);
    }

    /// Returns an immutable reference to the underlying bit data
    #[inline]
    pub fn as_slice(&self) -> &crate::data::BitSlice {
        self.0.row(0)
    }

    /// Returns a mutable reference to the underlying bit data
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut crate::data::BitSlice {
        self.0.row_mut(0)
    }

    /// Returns a reference to the underlying `BitMatrix`
    #[inline]
    pub fn as_matrix(&self) -> &BitMatrix {
        &self.0
    }

    /// Alias for `as_matrix`, emphasizing that the `BitVector` is being treated as a
    /// row vector, e.g. for matrix multiplication
    #[inline]
    pub fn as_row_vector(&self) -> &BitMatrix {
        &self.0
    }
}

/// Formats the vector for display
impl fmt::Display for BitVector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[")?;
        for i in 0..self.len() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", if self[i] { 1 } else { 0 })?;
        }
        write!(f, "]")
    }
}

/// XOR operation for BitVector
impl BitXor for &BitVector {
    type Output = BitVector;

    fn bitxor(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.len(),
            rhs.len(),
            "BitVectors must have the same length for XOR"
        );
        let mut result = self.clone();
        result.xor_with(rhs);
        result
    }
}

/// XOR operation for owned BitVector
impl BitXor for BitVector {
    type Output = BitVector;

    fn bitxor(mut self, rhs: Self) -> Self::Output {
        self ^= rhs;
        self
    }
}

/// XOR-assign operation for BitVector
impl BitXorAssign<&BitVector> for BitVector {
    fn bitxor_assign(&mut self, rhs: &BitVector) {
        self.xor_with(rhs);
    }
}

/// XOR-assign operation for owned BitVector
impl BitXorAssign<BitVector> for BitVector {
    fn bitxor_assign(&mut self, rhs: BitVector) {
        self.xor_with(&rhs);
    }
}

/// Equality comparison for BitVector
impl PartialEq for BitVector {
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            return false;
        }
        for i in 0..self.len() {
            if self.bit(i) != other.bit(i) {
                return false;
            }
        }
        true
    }
}

impl Eq for BitVector {}

/// Allows indexing into the vector to return the bit at `index`
///
/// `matrix[i]` is equivalent to `matrix.bit(i)`. Note this differs from how indexing works
/// for [`BitVec`], which indexes over [`BitBlock`]s, not individual bits.
impl Index<usize> for BitVector {
    type Output = bool;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        if self.bit(index) {
            &true
        } else {
            &false
        }
    }
}

impl From<BitVector> for BitMatrix {
    fn from(vector: BitVector) -> Self {
        vector.0
    }
}

impl TryFrom<BitMatrix> for BitVector {
    type Error = &'static str;

    fn try_from(matrix: BitMatrix) -> Result<Self, Self::Error> {
        if matrix.rows() != 1 {
            return Err("Cannot convert BitMatrix to BitVector unless it has exactly one row");
        }
        Ok(BitVector(matrix))
    }
}

impl Mul<&BitVector> for &BitMatrix {
    type Output = BitVector;

    fn mul(self, rhs: &BitVector) -> Self::Output {
        self.try_mul_vector(rhs).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::SmallRng, SeedableRng};

    #[test]
    fn test_matrix_vector_multiplication() {
        // Test basic matrix-vector multiplication
        // [1 0 1]   [1]   [1]
        // [0 1 1] * [1] = [1]
        // [1 1 0]   [0]   [0]
        //
        // Row 0: 1*1 + 0*1 + 1*0 = 1 + 0 + 0 = 1 (in GF(2))
        // Row 1: 0*1 + 1*1 + 1*0 = 0 + 1 + 0 = 1 (in GF(2))
        // Row 2: 1*1 + 1*1 + 0*0 = 1 + 1 + 0 = 0 (in GF(2), since 1⊕1=0)
        let matrix = BitMatrix::from_bool_vec(&vec![
            vec![true, false, true],
            vec![false, true, true],
            vec![true, true, false],
        ]);
        let vector = BitVector::from_bool_vec(&vec![true, true, false]);

        let result = &matrix * &vector;
        assert_eq!(result.len(), 3);
        assert!(result[0]); // 1*1 ⊕ 0*1 ⊕ 1*0 = 1
        assert!(result[1]); // 0*1 ⊕ 1*1 ⊕ 1*0 = 1
        assert!(!result[2]); // 1*1 ⊕ 1*1 ⊕ 0*0 = 0
    }

    #[test]
    fn test_matrix_vector_multiplication_identity() {
        // Test multiplication with identity matrix
        let identity = BitMatrix::identity(3);
        let vector = BitVector::from_bool_vec(&vec![true, false, true]);

        let result = &identity * &vector;
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], true);
        assert_eq!(result[1], false);
        assert_eq!(result[2], true);
    }

    #[test]
    fn test_matrix_vector_multiplication_zeros() {
        // Test multiplication with zero matrix
        let zero_matrix = BitMatrix::zeros(2, 3);
        let vector = BitVector::from_bool_vec(&vec![true, true, true]);

        let result = &zero_matrix * &vector;
        assert_eq!(result.len(), 2);
        assert!(result.is_zero());
    }

    #[test]
    #[should_panic(expected = "Cannot multiply matrix")]
    fn test_matrix_vector_multiplication_dimension_mismatch() {
        let matrix = BitMatrix::zeros(2, 3);
        let vector = BitVector::zeros(2); // Wrong dimension

        let _result = &matrix * &vector;
    }

    #[test]
    fn test_matrix_vector_multiplication_random() {
        let mut rng = SmallRng::seed_from_u64(123);
        let matrix = BitMatrix::random(&mut rng, 5, 7);
        let vector = BitVector::random(&mut rng, 7);

        let result = &matrix * &vector;
        assert_eq!(result.len(), 5);

        // Verify the multiplication manually for first row
        let mut expected_first = false;
        for j in 0..7 {
            expected_first ^= matrix.bit(0, j) & vector[j];
        }
        assert_eq!(result[0], expected_first);
    }
}
