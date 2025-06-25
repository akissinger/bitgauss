use std::{fmt, ops::Div};

use crate::{BitBlock, BitMatrix, BitVec};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ECCError(pub String);

// Standard implementations of error traits for `BitMatrixError`
impl std::fmt::Display for ECCError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ECCError: {}", self.0)
    }
}
impl std::error::Error for ECCError {}

pub struct BinaryLinearCode {
    pub(crate) n_codeword_length: usize,
    pub(crate) k_codespace_dimension: usize,
    pub(crate) d_code_distance: usize,
    pub(crate) g_generator_matrix: BitMatrix,
    known_is_standard: Option<bool>,
}

impl BinaryLinearCode {
    /// Each row of `to_encode` represents a bitvector being encoded
    /// That is to say `to_encode` is R by `k_codespace_dimension`
    /// and after encoding we get a result of R by `n_codeword_length`
    ///
    /// # Errors
    ///
    /// If the number of columns of `to_encode` are incorrect for the dimension
    /// of the codespace
    pub fn encode(&self, to_encode: &BitMatrix) -> Result<BitMatrix, ECCError> {
        if self.k_codespace_dimension != to_encode.cols() {
            return Err(ECCError(
                "the quantity to be encoded did not have the right number of columns".to_owned(),
            ));
        }
        Ok(to_encode * &self.g_generator_matrix)
    }

    /// the standard `[n,k,d]_q` notation
    pub fn n_k_d_q(&self) -> [usize; 4] {
        [
            self.n_codeword_length,
            self.k_codespace_dimension,
            self.d_code_distance,
            2,
        ]
    }

    /// Give the rate of information transfer.
    /// For a codespace of dimension `k` and codewords of length `n`
    /// the information transfer is slowed by a factor of `k/n`.
    /// This is compensated by the benefit of being able to detect and correct errors.
    pub fn rate<T: From<usize> + Div<T, Output = T>>(&self) -> T {
        let [n, k, _, _] = self.n_k_d_q();
        let n_t: T = n.into();
        let k_t: T = k.into();
        k_t / n_t
    }

    /// Is the G matrix in standard form with [I | P]
    ///
    /// # Panics
    ///
    /// The `k_codespace_dimension` by `n_codeword_length` matrix `g_generator_matrix`
    /// is assumed to have no more rows than columns
    /// the rate will be `<= 1`
    pub fn is_standard(&self) -> bool {
        if let Some(known_is_standard) = self.known_is_standard {
            return known_is_standard;
        }
        let cols_selected = (0..self.k_codespace_dimension).collect::<Vec<_>>();
        let left_block = self
            .g_generator_matrix
            .select_cols(&cols_selected)
            .expect("The matrix has fewer rows than columns by assumption");
        left_block == BitMatrix::identity(self.k_codespace_dimension)
    }
}

#[repr(transparent)]
pub struct HadamardCode(BinaryLinearCode);

impl HadamardCode {
    /// Each row of `to_encode` represents a bitvector being encoded
    /// That is to say `to_encode` is R by `k_codespace_dimension`
    /// and after encoding we get a result of R by `n_codeword_length`
    ///
    /// # Errors
    ///
    /// If the number of columns of `to_encode` are incorrect for the dimension
    /// of the codespace
    pub fn encode(&self, to_encode: &BitMatrix) -> Result<BitMatrix, ECCError> {
        self.0.encode(to_encode)
    }

    /// the standard `[n,k,d]_q` notation
    pub fn n_k_d_q(&self) -> [usize; 4] {
        self.0.n_k_d_q()
    }

    /// Give the rate of information transfer.
    /// For a codespace of dimension `k` and codewords of length `n`
    /// the information transfer is slowed by a factor of `k/n`.
    /// This is compensated by the benefit of being able to detect and correct errors.
    pub fn rate<T: From<usize> + Div<T, Output = T>>(&self) -> T {
        self.0.rate()
    }
}

impl HadamardCode {
    pub fn new(k_codespace_dimension: usize, augmented: bool) -> Self {
        let n_codeword_length: usize = 1 << k_codespace_dimension;
        let starting: BitBlock = if augmented {
            (n_codeword_length >> 1) as BitBlock
        } else {
            0
        };
        let iter = (starting..n_codeword_length as BitBlock).map(|z| {
            let b: BitVec = vec![z.rotate_right(k_codespace_dimension as u32)].into();
            BitMatrix::row_vector(k_codespace_dimension, b)
        });
        let mut g_generator_matrix = BitMatrix::vstack_from_owned_iter(iter);
        g_generator_matrix.transpose_inplace();
        let n_codeword_length = if augmented {
            n_codeword_length >> 1
        } else {
            n_codeword_length
        };
        Self(BinaryLinearCode {
            n_codeword_length,
            k_codespace_dimension,
            d_code_distance: n_codeword_length >> 1,
            g_generator_matrix,
            known_is_standard: Some(false),
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn hadamard_3_8() {
        let cur_had = HadamardCode::new(3, false);
        assert_eq!(cur_had.0.n_codeword_length, 8);
        assert_eq!(cur_had.0.k_codespace_dimension, 3);
        assert_eq!(cur_had.0.d_code_distance, 4);
        assert_eq!(cur_had.0.g_generator_matrix.rows(), 3);
        assert_eq!(cur_had.0.g_generator_matrix.cols(), 8);
        let a = [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1],
        ];
        for k_idx in 0..3 {
            for n_idx in 0..8 {
                assert_eq!(
                    a[n_idx][k_idx] == 1,
                    cur_had.0.g_generator_matrix.bit(k_idx, n_idx),
                    "G[{k_idx}{n_idx}]"
                );
            }
        }
        assert_eq!(
            cur_had.0.g_generator_matrix,
            BitMatrix::build(3, 8, |i, j| a[j][i] == 1)
        );
    }

    #[test]
    fn augmented_hadamard_3_4() {
        let cur_had = HadamardCode::new(3, true);
        assert_eq!(cur_had.0.n_codeword_length, 4);
        assert_eq!(cur_had.0.k_codespace_dimension, 3);
        assert_eq!(cur_had.0.d_code_distance, 2);
        assert_eq!(cur_had.0.g_generator_matrix.rows(), 3);
        assert_eq!(cur_had.0.g_generator_matrix.cols(), 4);
        let a = [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]];
        for k_idx in 0..3 {
            for n_idx in 0..4 {
                assert_eq!(
                    a[n_idx][k_idx] == 1,
                    cur_had.0.g_generator_matrix.bit(k_idx, n_idx),
                    "G[{k_idx}{n_idx}]"
                );
            }
        }
        assert_eq!(
            cur_had.0.g_generator_matrix,
            BitMatrix::build(3, 4, |i, j| a[j][i] == 1)
        );
    }
}
