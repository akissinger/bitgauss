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

enum ECCPresentations {
    JustG(BitMatrix),
    StandardP(BitMatrix),
    #[allow(dead_code)]
    BothGH(BitMatrix, BitMatrix),
}

impl ECCPresentations {
    #[allow(dead_code)]
    fn valid(&self) -> bool {
        match self {
            ECCPresentations::JustG(bit_matrix) => bit_matrix.rows() < bit_matrix.cols(),
            ECCPresentations::StandardP(_) => true,
            ECCPresentations::BothGH(g_matrix, h_matrix) => {
                let (k_g, n_g) = (g_matrix.rows(), g_matrix.cols());
                let (n_minus_k_h, n_h) = (h_matrix.rows(), h_matrix.cols());
                if n_g == n_h && n_minus_k_h + k_g == n_g {
                    let g_transpose = g_matrix.transposed();
                    let h_times_g_t = h_matrix * &g_transpose;
                    h_times_g_t == BitMatrix::zeros(n_minus_k_h, k_g)
                } else {
                    false
                }
            }
        }
    }
}

pub struct BinaryLinearCode {
    n_codeword_length: usize,
    k_codespace_dimension: usize,
    d_code_distance: usize,
    chosen_presentation: ECCPresentations,
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
        match &self.chosen_presentation {
            ECCPresentations::JustG(g_generator_matrix)
            | ECCPresentations::BothGH(g_generator_matrix, _) => Ok(to_encode * g_generator_matrix),
            ECCPresentations::StandardP(p_part_of_g) => {
                let left_part = to_encode.clone();
                let right_part = to_encode * p_part_of_g;
                Ok(BitMatrix::hstack(&left_part, &right_part))
            }
        }
    }

    /// Each column of `to_decode` represents a bitvector being decoded
    /// That is to say `to_decode` is `n_codeword_length` by C
    /// and after decoding we would get a result of `k_codespace_dimension` by C
    ///
    /// # Errors
    ///
    /// If the number of rows of `to_decode` are incorrect for the dimension
    /// of the codewords
    pub fn is_codeword(&self, to_decode: &BitMatrix) -> Result<BitVec, ECCError> {
        if self.n_codeword_length != to_decode.rows() {
            return Err(ECCError(
                "the quantity to be decoded did not have the right number of rows".to_owned(),
            ));
        }
        match &self.chosen_presentation {
            ECCPresentations::JustG(_) => Err(ECCError(
                "We do not have the check matrix, we need to compute it".to_owned(),
            )),
            ECCPresentations::StandardP(bit_matrix) => {
                let p_transpose = bit_matrix.transposed();
                let identity_part = BitMatrix::identity(p_transpose.rows());
                let full_h = BitMatrix::hstack(&p_transpose, &identity_part);
                let count_columns = to_decode.cols();
                let mut to_check_zero_columns = &full_h * to_decode;
                to_check_zero_columns.transpose_inplace();
                let to_check_zero_rows = to_check_zero_columns;
                Ok((0..count_columns)
                    .map(|cur_row| to_check_zero_rows.row(cur_row).has_any_ones())
                    .collect())
            }
            ECCPresentations::BothGH(_, full_h) => {
                let count_columns = to_decode.cols();
                let mut to_check_zero_columns = full_h * to_decode;
                to_check_zero_columns.transpose_inplace();
                let to_check_zero_rows = to_check_zero_columns;
                Ok((0..count_columns)
                    .map(|cur_row| to_check_zero_rows.row(cur_row).has_any_ones())
                    .collect())
            }
        }
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

    /// Is the G matrix in standard form with [I | P]?
    /// This can be either explicitly as the `StandardP` variant or possibly checking the left block
    /// to see if it is the identity matrix.
    ///
    /// # Panics
    ///
    /// The `k_codespace_dimension` by `n_codeword_length` matrix `g_generator_matrix`
    /// is assumed to have no more rows than columns
    /// the rate will be `<= 1`
    pub fn is_standard(&self, only_obviously: bool) -> bool {
        if only_obviously {
            return matches!(self.chosen_presentation, ECCPresentations::StandardP(_));
        }
        match &self.chosen_presentation {
            ECCPresentations::JustG(g_generator_matrix)
            | ECCPresentations::BothGH(g_generator_matrix, _) => {
                let cols_selected = (0..self.k_codespace_dimension).collect::<Vec<_>>();
                let left_block = g_generator_matrix
                    .select_cols(&cols_selected)
                    .expect("The matrix has fewer rows than columns by assumption");
                left_block == BitMatrix::identity(self.k_codespace_dimension)
            }
            ECCPresentations::StandardP(_) => true,
        }
    }

    /// Either give a reference to the `G` matrix in the `Ok` variant
    /// or reconstruct it from the `P` matrix in the `Err` variant
    #[allow(clippy::missing_errors_doc)]
    pub fn g_generator_matrix(&self) -> Result<&BitMatrix, BitMatrix> {
        match &self.chosen_presentation {
            ECCPresentations::JustG(bit_matrix) | ECCPresentations::BothGH(bit_matrix, _) => {
                Ok(bit_matrix)
            }
            ECCPresentations::StandardP(bit_matrix) => Err(BitMatrix::hstack(
                &BitMatrix::identity(self.k_codespace_dimension),
                bit_matrix,
            )),
        }
    }

    #[allow(clippy::missing_panics_doc)]
    pub fn standardize(self) -> Self {
        if self.is_standard(true) {
            return self;
        }
        match self.chosen_presentation {
            ECCPresentations::JustG(mut bit_matrix)
            | ECCPresentations::BothGH(mut bit_matrix, _) => {
                let pivot_cols = bit_matrix.gauss(true);
                let num_cols = bit_matrix.cols();
                let chosen_presentation = if pivot_cols.len() == self.k_codespace_dimension {
                    let cols_for_p_selected = (0..num_cols)
                        .filter(|z| !pivot_cols.contains(z))
                        .collect::<Vec<_>>();
                    let p_matrix = bit_matrix
                        .select_cols(&cols_for_p_selected)
                        .expect("Manifestly all within range of logical columns");
                    ECCPresentations::StandardP(p_matrix)
                } else {
                    ECCPresentations::JustG(bit_matrix)
                };
                Self {
                    n_codeword_length: self.n_codeword_length,
                    k_codespace_dimension: self.k_codespace_dimension,
                    d_code_distance: self.d_code_distance,
                    chosen_presentation,
                }
            }
            ECCPresentations::StandardP(_) => {
                unreachable!("Already checked that it was not standard already");
            }
        }
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
            chosen_presentation: ECCPresentations::JustG(g_generator_matrix),
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
        assert_eq!(
            cur_had
                .0
                .g_generator_matrix()
                .expect("Is given as G")
                .rows(),
            3
        );
        assert_eq!(
            cur_had
                .0
                .g_generator_matrix()
                .expect("Is given as G")
                .cols(),
            8
        );
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
                    cur_had
                        .0
                        .g_generator_matrix()
                        .expect("Is given as G")
                        .bit(k_idx, n_idx),
                    "G[{k_idx}{n_idx}]"
                );
            }
        }
        assert_eq!(
            *cur_had.0.g_generator_matrix().expect("Is given as G"),
            BitMatrix::build(3, 8, |i, j| a[j][i] == 1)
        );
    }

    #[test]
    fn augmented_hadamard_3_4() {
        let cur_had = HadamardCode::new(3, true);
        assert_eq!(cur_had.0.n_codeword_length, 4);
        assert_eq!(cur_had.0.k_codespace_dimension, 3);
        assert_eq!(cur_had.0.d_code_distance, 2);
        assert_eq!(
            cur_had
                .0
                .g_generator_matrix()
                .expect("Is given as G")
                .rows(),
            3
        );
        assert_eq!(
            cur_had
                .0
                .g_generator_matrix()
                .expect("Is given as G")
                .cols(),
            4
        );
        let a = [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]];
        for k_idx in 0..3 {
            for n_idx in 0..4 {
                assert_eq!(
                    a[n_idx][k_idx] == 1,
                    cur_had
                        .0
                        .g_generator_matrix()
                        .expect("Is given as G")
                        .bit(k_idx, n_idx),
                    "G[{k_idx}{n_idx}]"
                );
            }
        }
        assert_eq!(
            *cur_had.0.g_generator_matrix().expect("Is given as G"),
            BitMatrix::build(3, 4, |i, j| a[j][i] == 1)
        );
    }
}
