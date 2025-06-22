use crate::bitvec::*;
use rand::Rng;
use std::{
    fmt,
    ops::{Index, Mul},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitMatrixError(pub String);

// Standard implementations of error traits for `BitMatrixError`
impl std::fmt::Display for BitMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BitMatrixError: {}", self.0)
    }
}
impl std::error::Error for BitMatrixError {}

/// A matrix of bits, represented as a vector of blocks of bits
///
/// The matrix is stored in row-major order, with each row represented as a `BitRange` of `BitBlock`s. If
/// the number of columns is not a multiple of [`BLOCKSIZE`], the last block in each row will be padded with 0s.
///
/// The matrix is additionally allowed to be padded arbitrarily in either dimension, e.g. to make it square. In
/// that case, extra bits beyond `rows` and `cols` should always be 0.
///
/// The logical rows and columns of the matrix are given by `rows` and `cols`, while the full padded number of
/// rows is given by `data.len() / col_blocks` and the padded number of columns is `col_blocks * BLOCKSIZE`.
#[derive(Clone, Debug)]
pub struct BitMatrix {
    /// the number of logical rows in the matrix
    rows: usize,

    /// the number of logical columns in the matrix
    cols: usize,

    /// the number of [`BitBlock`]s used to store each row, i.e. actual 2D matrix has `col_blocks * BLOCKSIZE` many columns
    col_blocks: usize,

    /// a [`BitVec`] containing the data of the matrix, stored in row-major order
    data: BitVec,
}

/// A trait for types that can have row operations performed on them
pub trait RowOps {
    fn add_row(&mut self, from: usize, to: usize);
    fn swap_rows(&mut self, from: usize, to: usize);
}

impl BitMatrix {
    /// Gets the bit at position `(i, j)`
    #[inline]
    pub fn bit(&self, i: usize, j: usize) -> bool {
        self.data.bit(self.col_blocks * BLOCKSIZE * i + j)
    }

    /// Sets the bit at position `(i, j)` to `b`
    #[inline]
    pub fn set_bit(&mut self, i: usize, j: usize, b: bool) {
        self.data.set_bit(self.col_blocks * BLOCKSIZE * i + j, b);
    }

    /// Builds a `BitMatrix` from a function `f` that determines the value of each bit
    ///
    /// # Arguments
    /// * `rows` - the number of rows in the matrix
    /// * `cols` - the number of columns in the matrix
    /// * `f` - a function that takes the row and column indices and returns a boolean value for each bit
    pub fn build(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> bool) -> Self {
        let col_blocks = min_blocks(cols);
        let data = (0..rows)
            .flat_map(|i| (0..BLOCKSIZE * col_blocks).map(move |j| (i, j)))
            .map(|(i, j)| if i < rows && j < cols { f(i, j) } else { false })
            .collect();
        BitMatrix {
            rows,
            cols,
            col_blocks,
            data,
        }
    }

    /// Creates a new `BitMatrix` of size `rows` x `cols` with all bits set to 0
    pub fn zeros(rows: usize, cols: usize) -> Self {
        let col_blocks = min_blocks(cols);
        BitMatrix {
            rows,
            cols,
            col_blocks,
            data: BitVec::zeros(rows * col_blocks),
        }
    }

    /// Checks if the matrix consists of all zero bits
    pub fn is_zero(&self) -> bool {
        self.data.is_empty() || self.data.block_iter().all(|b| b == 0)
    }

    /// Creates a new identity `BitMatrix` of size `size` x `size`
    pub fn identity(size: usize) -> Self {
        let blocks = min_blocks(size);
        let num_blocks = size * blocks;

        let data = (0..num_blocks)
            .map(|i| {
                let row = i / blocks;
                let col_block = i % blocks;
                if row / BLOCKSIZE == col_block && i < size * size {
                    MSB_ON >> (row % BLOCKSIZE)
                } else {
                    0
                }
            })
            .collect();
        BitMatrix {
            rows: size,
            cols: size,
            col_blocks: blocks,
            data,
        }
    }

    /// Creates a new random `BitMatrix` of size `rows` x `cols`
    ///
    /// Bits outside of the logical size of the matrix (i.e. `rows` and `cols`) will be masked to 0.
    #[inline]
    pub fn random(rng: &mut impl Rng, rows: usize, cols: usize) -> Self {
        let col_blocks = min_blocks(cols);
        let num_blocks = rows * col_blocks;
        let mask = BitBlock::MAX.wrapping_shl((BLOCKSIZE - (cols % BLOCKSIZE)) as u32);
        let data = (0..num_blocks)
            .map(|i| {
                if i % col_blocks == col_blocks - 1 {
                    mask & rng.random::<BitBlock>()
                } else {
                    rng.random::<BitBlock>()
                }
            })
            .collect();
        BitMatrix {
            rows,
            cols,
            col_blocks,
            data,
        }
    }

    /// Creates a new random invertible `BitMatrix` of size `size` x `size`
    #[inline]
    pub fn random_invertible(rng: &mut impl Rng, size: usize) -> Self {
        let mut m = BitMatrix::identity(size);

        for _ in 0..10 * size * size {
            let r1 = rng.random_range(0..size);
            let mut r2 = rng.random_range(0..size - 1);
            if r2 >= r1 {
                r2 += 1;
            }
            m.add_row(r1, r2);
        }

        m
    }

    /// Returns the number of logical rows in the matrix
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Returns the number of logical columns in the matrix
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Adds (XORs) the bits from a `BitRange` to a specified row
    #[inline]
    pub fn add_bits_to_row(&mut self, bits: &BitSlice, row: usize) {
        self.data.xor_in(bits, row * self.col_blocks);
    }

    /// Returns an immutable reference to a row of the matrix as a `BitRange`
    #[inline]
    pub fn row(&self, row: usize) -> &BitSlice {
        &self.data[row * self.col_blocks..(row + 1) * self.col_blocks]
    }

    /// Returns a mutable reference to a row of the matrix as a `BitRange`
    #[inline]
    pub fn row_mut(&mut self, row: usize) -> &mut BitSlice {
        &mut self.data[row * self.col_blocks..(row + 1) * self.col_blocks]
    }

    /// Pads the matrix with zero-rows and zero-columns to make it square
    #[inline]
    pub fn pad_to_square(&mut self) {
        let data_rows = self.data.len() / self.col_blocks;
        let row_blocks = min_blocks(data_rows);
        if data_rows != row_blocks * BLOCKSIZE || row_blocks != self.col_blocks {
            let blocks = usize::max(row_blocks, self.col_blocks);
            let mut data = BitVec::with_capacity(BLOCKSIZE * blocks * blocks);
            for i in 0..(BLOCKSIZE * blocks) {
                for j in 0..blocks {
                    data.push_block(if i < self.rows() && j < self.col_blocks {
                        self.data[i * self.col_blocks + j]
                    } else {
                        0
                    });
                }
            }

            self.data = data;
            self.col_blocks = blocks;
        }
    }

    /// Main working function for transposition
    ///
    /// If `source` is given as `Some(bit_matrix)`, then copy bits from `bit_matrix` in transposed
    /// position. Otherwise if it is `None` then transpose in place (assuming the matrix is already
    /// padded to be square).
    fn transpose_helper(&mut self, source: Option<&BitMatrix>) {
        let mut buffer: [BitBlock; BLOCKSIZE] = [0; BLOCKSIZE];
        for i in 0..min_blocks(self.rows) {
            for j in 0..self.col_blocks {
                let dest_block = BLOCKSIZE * i * self.col_blocks + j;
                let source_block;
                if let Some(m) = source {
                    source_block = BLOCKSIZE * j * m.col_blocks + i;
                    // load source_block into buffer
                    for k in 0..BLOCKSIZE {
                        let l = source_block + k * m.col_blocks;
                        buffer[k] = if l < m.data.len() { m.data[l] } else { 0 };
                    }
                } else {
                    source_block = BLOCKSIZE * j * self.col_blocks + i;
                    for k in 0..BLOCKSIZE {
                        // if this block is above the diagonal, swap it with the one in transposed position
                        if i < j {
                            self.data.swap(
                                source_block + k * self.col_blocks,
                                dest_block + k * self.col_blocks,
                            );
                        }

                        // load dest_block into buffer
                        buffer[k] = self.data[dest_block + k * self.col_blocks];
                    }
                }

                // transpose the block in place by iteratively transposing blocks of half the size
                // until we get down to block size 1
                let mut swap_width = BLOCKSIZE;
                let mut swap_mask0 = BitBlock::MAX;
                while swap_width != 1 {
                    swap_width >>= 1;

                    // masks that pick the left half of the bits and right half of the bits in each block
                    swap_mask0 ^= swap_mask0 >> swap_width;
                    let swap_mask1 = BitBlock::MAX ^ swap_mask0;

                    for block_row in (0..BLOCKSIZE).step_by(swap_width * 2) {
                        for row in block_row..block_row + swap_width {
                            let b0 = buffer[row];
                            let b1 = buffer[row + swap_width];
                            buffer[row] = (b0 & swap_mask0) | ((b1 & swap_mask0) >> swap_width);
                            buffer[row + swap_width] =
                                (b1 & swap_mask1) | ((b0 & swap_mask1) << swap_width);
                        }
                    }
                }

                for k in 0..BLOCKSIZE {
                    let l = dest_block + k * self.col_blocks;
                    if l < self.data.len() {
                        self.data[l] = buffer[k];
                    }
                }
            }
        }
    }

    /// Returns a transposed copy of the matrix
    #[inline]
    pub fn transposed(&self) -> Self {
        let mut dest = Self::zeros(self.cols, self.rows);
        dest.transpose_helper(Some(self));
        dest
    }

    /// Transposes the matrix in place, padding allocated memory with 0s if necessary
    #[inline]
    pub fn transpose_inplace(&mut self) {
        self.pad_to_square();
        (self.rows, self.cols) = (self.cols, self.rows);
        self.transpose_helper(None);
    }

    /// Performs gaussian elimination while also performing matching row operations on `proxy`
    /// and returns a vector of pivot columns.
    fn gauss_helper(&mut self, full: bool, proxy: &mut impl RowOps) -> Vec<usize> {
        let mut row = 0;
        let mut pcol = 0;
        let mut pcols = vec![];
        while row < self.rows() {
            let mut next_row = None;
            'outer: while pcol < self.cols() {
                for i in row..self.rows() {
                    if self[(i, pcol)] {
                        next_row = Some(i);
                        break 'outer;
                    }
                }
                pcol += 1;
            }

            if let Some(row1) = next_row {
                if row != row1 {
                    self.swap_rows(row, row1);
                    proxy.swap_rows(row, row1);
                }

                let row_vec = self.row(row).to_vec();

                for i in (row1 + 1)..self.rows() {
                    if self[(i, pcol)] {
                        self.add_bits_to_row(&row_vec, i);
                        proxy.add_row(row, i);
                    }
                }

                row += 1;
                pcols.push(pcol);
                pcol += 1;
            } else {
                break;
            }
        }

        if full {
            for row in (0..pcols.len()).rev() {
                let pcol = pcols[row];
                let row_vec = self.row(row).to_vec();
                for i in 0..row {
                    if self[(i, pcol)] {
                        self.add_bits_to_row(&row_vec, i);
                        proxy.add_row(row, i);
                    }
                }
            }
        }

        pcols
    }

    /// Performs gaussian elimination
    ///
    /// If `full` is true, then perform full Gauss-Jordan to produce reduced echelon form, otherwise
    /// just returns echelon form.
    #[inline]
    pub fn gauss(&mut self, full: bool) {
        self.gauss_helper(full, &mut ());
    }

    /// Computes the rank of the matrix using gaussian elimination
    #[inline]
    pub fn rank(&self) -> usize {
        self.clone().gauss_helper(false, &mut ()).len()
    }

    /// Computes the inverse of the matrix if it is invertible, otherwise returns an error
    pub fn try_inverse(&self) -> Result<Self, BitMatrixError> {
        if self.rows() != self.cols() {
            return Err(BitMatrixError("Matrix must be square".to_string()));
        }
        let mut inv = BitMatrix::identity(self.cols());
        let pcols = self.clone().gauss_helper(true, &mut inv);

        if pcols.len() != self.cols() {
            return Err(BitMatrixError("Matrix is not invertible".to_string()));
        }

        Ok(inv)
    }

    /// Computes the inverse of an invertible matrix
    pub fn inverse(&self) -> Self {
        self.try_inverse().unwrap()
    }

    /// Tries to multiply two matrices and returns the result
    ///
    /// Returns an error if the matrices have incompatible dimensions
    pub fn try_mul(&self, other: &Self) -> Result<Self, BitMatrixError> {
        if self.cols() != other.rows() {
            return Err(BitMatrixError(format!(
                "Cannot multiply matrices of dimensions {}x{} and {}x{}",
                self.rows(),
                self.cols(),
                other.rows(),
                other.cols()
            )));
        }

        let mut res = BitMatrix::zeros(self.rows, other.cols);

        for i in 0..self.rows {
            let row = res.row_mut(i);
            self.row(i).iter().enumerate().for_each(|(j, b)| {
                if b {
                    *row ^= other.row(j);
                }
            });
        }

        Ok(res)
    }

    /// Try to vertically stack this matrix with another one and returns the result
    ///
    /// The resulting matrix will have the minimal column padding. Returns an error if the
    /// matrices have different numbers of columns.
    pub fn try_vstack(&self, other: &Self) -> Result<Self, BitMatrixError> {
        if self.cols() != other.cols() {
            return Err(BitMatrixError(format!(
                "Cannot vertically stack matrices with different number of columns: {} != {}",
                self.cols(),
                other.cols()
            )));
        }

        let rows = self.rows() + other.rows();
        let mut data = BitVec::with_capacity(rows * self.col_blocks);
        let col_blocks = min_blocks(self.cols());

        data.reserve(rows * col_blocks);
        for i in 0..self.rows() {
            let start = i * self.col_blocks;
            data.extend_from_slice(&self.data[start..start + col_blocks]);
        }

        for i in 0..other.rows() {
            let start = i * other.col_blocks;
            data.extend_from_slice(&other.data[start..start + col_blocks]);
        }

        Ok(BitMatrix {
            rows,
            cols: self.cols(),
            col_blocks,
            data,
        })
    }

    /// Vertically stacks this matrix with another one and returns the result
    pub fn vstack(&self, other: &Self) -> Self {
        self.try_vstack(other).unwrap()
    }

    /// Horizontally stacks this matrix with another one and returns the result
    pub fn try_hstack(&self, other: &Self) -> Result<Self, BitMatrixError> {
        if self.rows() != other.rows() {
            return Err(BitMatrixError(format!(
                "Cannot horizontally stack matrices with different number of rows: {} != {}",
                self.rows(),
                other.rows()
            )));
        }

        let cols = self.cols() + other.cols();
        let mut data = BitVec::with_capacity(self.rows * min_blocks(cols));
        let col_blocks = min_blocks(cols);
        let self_col_blocks = min_blocks(self.cols());
        let other_col_blocks = min_blocks(other.cols());
        let pop_one = self_col_blocks + other_col_blocks > col_blocks;
        let shift = BLOCKSIZE * self_col_blocks - self.cols();

        for i in 0..self.rows() {
            let start_self = i * self.col_blocks;
            let start_other = i * other.col_blocks;
            data.extend_from_slice(&self.data[start_self..start_self + self_col_blocks]);
            data.extend_from_slice_left_shifted(
                &other.data[start_other..start_other + other_col_blocks],
                shift,
            );

            // pop the last block from the row if the extra padding is not needed
            if pop_one {
                data.pop();
            }
        }

        Ok(BitMatrix {
            rows: self.rows(),
            cols,
            col_blocks,
            data,
        })
    }

    /// Horizontally stacks this matrix with another one and returns the result
    pub fn hstack(&self, other: &Self) -> Self {
        self.try_hstack(other).unwrap()
    }

    /// Vertically stacks an iterator of `BitMatrix` instances into a single `BitMatrix`
    ///
    /// If the iterator is empty, returns an empty `BitMatrix` with 0 rows and 0 columns.
    pub fn vstack_from_iter<'a>(iter: impl IntoIterator<Item = &'a BitMatrix>) -> Self {
        let mut it = iter.into_iter();
        if let Some(first) = it.next() {
            it.fold(first.clone(), |m, n| m.vstack(n))
        } else {
            Self::zeros(0, 0)
        }
    }

    /// Horizontally stacks an iterator of `BitMatrix` instances into a single `BitMatrix`
    ///
    /// If the iterator is empty, returns an empty `BitMatrix` with 0 rows and 0 columns.
    pub fn hstack_from_iter<'a>(iter: impl IntoIterator<Item = &'a BitMatrix>) -> Self {
        let mut it = iter.into_iter();
        if let Some(first) = it.next() {
            it.fold(first.clone(), |m, n| m.hstack(n))
        } else {
            Self::zeros(0, 0)
        }
    }

    /// Computes a basis for the nullspace a the matrix and returns it as the rows of a new matrix
    pub fn nullspace(&self) -> Vec<BitMatrix> {
        if self.rows() == 0 || self.cols() == 0 {
            return Vec::new();
        }

        let mut m = self.clone();
        let pivot_cols = m.gauss_helper(true, &mut ());
        let mut free_vars = Vec::with_capacity(self.cols() - pivot_cols.len());
        let mut it = pivot_cols.iter().peekable();

        for i in 0..self.cols() {
            if it.peek().is_some_and(|&&p| p == i) {
                it.next();
            } else {
                free_vars.push(i);
            }
        }

        // Generate basis vectors for the nullspace
        let mut basis = Vec::with_capacity(free_vars.len());
        for &free_var in &free_vars {
            let mut vec = Self::zeros(1, self.cols());
            vec.set_bit(0, free_var, true);

            // Back substitution
            for (row, &pivot_col) in pivot_cols.iter().enumerate().rev() {
                if free_var > pivot_col && m[(row, free_var)] {
                    vec.set_bit(0, pivot_col, true);
                }
            }

            basis.push(vec);
        }

        basis
    }
}

/// Two matrices are considered equal if they represent the same logical matrix, possibly with different
/// padding (i.e. col_blocks and row_blocks can be different)
impl PartialEq for BitMatrix {
    fn eq(&self, other: &Self) -> bool {
        if self.rows() != other.rows() || self.cols() != other.cols() {
            return false;
        }

        for i in 0..self.rows() {
            for j in 0..self.col_blocks {
                if j * BLOCKSIZE >= self.cols() {
                    break;
                } else if self.data[i * self.col_blocks + j] != other.data[i * other.col_blocks + j]
                {
                    return false;
                }
            }
        }

        true
    }
}

impl Eq for BitMatrix {}

/// A no-op implementation of `RowOps` for when we don't need to track row operations
impl RowOps for () {
    #[inline]
    fn add_row(&mut self, _: usize, _: usize) {}

    #[inline]
    fn swap_rows(&mut self, _: usize, _: usize) {}
}

/// An implementation of `RowOps` for `BitMatrix` that allows performing row operations on the matrix
impl RowOps for BitMatrix {
    #[inline]
    fn add_row(&mut self, from: usize, to: usize) {
        self.data.xor_range(
            from * self.col_blocks,
            to * self.col_blocks,
            self.col_blocks,
        );
    }

    #[inline]
    fn swap_rows(&mut self, from: usize, to: usize) {
        self.data.swap_range(
            from * self.col_blocks,
            to * self.col_blocks,
            self.col_blocks,
        );
    }
}

/// Allows indexing into the matrix to return the bit at `(row, col)
///
/// `matrix[(row, col)]` is equivalent to `matrix.bit(row, col)`. Note this differs from how indexing works for [`BitVec`], which indexes
/// over [`BitBlock`]s, not individual bits.
impl Index<(usize, usize)> for BitMatrix {
    type Output = bool;

    #[inline]
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        if self.bit(index.0, index.1) {
            &true
        } else {
            &false
        }
    }
}

/// Formats the matrix for display
///
/// Padding bits are not shown, only the logical size of the matrix is displayed.
impl fmt::Display for BitMatrix {
    /// Formats the matrix for display.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, " {} ", if self[(i, j)] { 1 } else { 0 })?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

impl Mul for &BitMatrix {
    type Output = BitMatrix;
    /// Multiplies two matrices.
    fn mul(self, rhs: Self) -> Self::Output {
        self.try_mul(rhs).unwrap()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{rngs::SmallRng, SeedableRng};

    #[test]
    fn random_gauss() {
        let mut rng = SmallRng::seed_from_u64(1);
        let mut m = BitMatrix::random(&mut rng, 8, 20); // will pad to 64x64

        println!(
            "rows {} cols {} blocks {}\n mask {:064b}",
            m.rows,
            m.cols,
            m.col_blocks,
            BitBlock::MAX.wrapping_shl((BLOCKSIZE - (m.cols % BLOCKSIZE)) as u32)
        );

        println!("{}", m);
        m.gauss(true);
        println!("{}", m);
    }

    #[test]
    fn identity() {
        let m = BitMatrix::identity(100);
        for i in 0..100 {
            for j in 0..100 {
                assert_eq!(m[(i, j)], i == j);
            }
        }
    }

    #[test]
    fn transpose() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 10, 4);
        let n = m.transposed();
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(j, i)]);
            }
        }

        let m = BitMatrix::random(&mut rng, 300, 200);
        let n = m.transposed();
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(j, i)]);
            }
        }
    }

    #[test]
    fn pad_to_square_sm() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 4, 5);
        let mut n = m.clone();
        n.pad_to_square();
        assert_eq!(n.col_blocks, 1);
        assert_eq!(n.data.len() / (n.col_blocks * BLOCKSIZE), 1);
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(i, j)]);
            }
        }
    }

    #[test]
    fn pad_to_square() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 300, 200);
        let mut n = m.clone();
        n.pad_to_square();
        assert_eq!(n.col_blocks, 5);
        assert_eq!(n.data.len() / (n.col_blocks * BLOCKSIZE), 5);
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(i, j)]);
            }
        }
    }

    #[test]
    fn transpose_inplace() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 10, 4);
        let mut n = m.clone();
        n.transpose_inplace();
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(j, i)]);
            }
        }
        n.transpose_inplace();
        assert_eq!(m, n);

        let m = BitMatrix::random(&mut rng, 300, 200);
        let mut n = m.clone();
        n.transpose_inplace();
        for i in 0..m.rows() {
            for j in 0..m.cols() {
                assert_eq!(m[(i, j)], n[(j, i)]);
            }
        }
        n.transpose_inplace();
        assert_eq!(m, n);
    }

    #[test]
    fn matrix_mult() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m1 = BitMatrix::random(&mut rng, 80, 100);
        let m2 = BitMatrix::random(&mut rng, 100, 70);
        let m3 = &m1 * &m2;

        for i in 0..m3.rows() {
            for j in 0..m3.cols() {
                let mut b = false;
                for k in 0..m1.cols() {
                    b ^= m1.bit(i, k) & m2.bit(k, j);
                }
                assert_eq!(m3.bit(i, j), b);
            }
        }
        // println!("{}\n*\n{}\n=\n{}", m1, m2, m3);
    }

    #[test]
    fn matrix_inv() {
        let mut rng = SmallRng::seed_from_u64(1);
        let sz = 100;
        let m = BitMatrix::random_invertible(&mut rng, sz);
        let n = m.inverse();
        let id = BitMatrix::identity(sz);

        assert_eq!(&m * &n, id);
        assert_eq!(&n * &m, id);
    }

    // test that two matrices with different padding but same logical size are equal
    #[test]
    fn matrix_eq_padding() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m1 = BitMatrix::random(&mut rng, 10, 20);
        let mut m2 = m1.clone();
        m2.pad_to_square();
        assert_eq!(m1, m2);
    }

    // test BitMatrix::vstack
    #[test]
    fn matrix_vstack() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m1 = BitMatrix::random(&mut rng, 10, 20);
        let m2 = BitMatrix::random(&mut rng, 5, 20);
        let m3 = m1.vstack(&m2);

        assert_eq!(m3.rows(), m1.rows() + m2.rows());
        assert_eq!(m3.cols(), m1.cols());
        for i in 0..m1.rows() {
            for j in 0..m1.cols() {
                assert_eq!(m3[(i, j)], m1[(i, j)]);
            }
        }
        for i in 0..m2.rows() {
            for j in 0..m2.cols() {
                assert_eq!(m3[(i + m1.rows(), j)], m2[(i, j)]);
            }
        }
    }

    // test BitMatrix::hstack
    #[test]
    fn matrix_hstack() {
        let cases = [(10, 20, 5), (10, 150, 5), (10, 200, 300)];

        for (rows, cols1, cols2) in cases {
            println!(
                "Testing hstack with {}x{} and {}x{}",
                rows, cols1, rows, cols2
            );
            let mut rng = SmallRng::seed_from_u64(1);
            let m1 = BitMatrix::random(&mut rng, rows, cols1);
            let m2 = BitMatrix::random(&mut rng, rows, cols2);
            let m3 = m1.hstack(&m2);

            assert_eq!(m3.rows(), m1.rows());
            assert_eq!(m3.cols(), m1.cols() + m2.cols());
            for i in 0..m1.rows() {
                for j in 0..m1.cols() {
                    assert_eq!(m3[(i, j)], m1[(i, j)]);
                }
                for j in 0..m2.cols() {
                    assert_eq!(m3[(i, j + m1.cols())], m2[(i, j)]);
                }
            }

            // check extra padding was not added
            assert_eq!(m3.col_blocks, min_blocks(m3.cols()));
        }
    }

    // test BitMatrix::nullspace
    #[test]
    fn matrix_nullspace() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 70, 200);
        let ns_mat = BitMatrix::vstack_from_iter(&m.nullspace());
        assert_eq!(ns_mat.rank(), ns_mat.rows());
        assert!((&m * &ns_mat.transposed()).is_zero());
    }

    #[test]
    fn test_build_function() {
        // Test building a matrix with a custom function
        let m = BitMatrix::build(3, 4, |i, j| (i + j) % 2 == 0);
        assert_eq!(m.rows(), 3);
        assert_eq!(m.cols(), 4);

        // Check the pattern: alternating bits
        for i in 0..3 {
            for j in 0..4 {
                assert_eq!(m[(i, j)], (i + j) % 2 == 0);
            }
        }
    }

    #[test]
    fn test_empty_matrix() {
        let m = BitMatrix::zeros(0, 0);
        assert_eq!(m.rows(), 0);
        assert_eq!(m.cols(), 0);
        assert!(m.is_zero());
    }

    #[test]
    fn test_single_element_matrix() {
        let m = BitMatrix::build(1, 1, |_, _| true);
        assert_eq!(m.rows(), 1);
        assert_eq!(m.cols(), 1);
        assert!(m[(0, 0)]);
        assert!(!m.is_zero());
    }

    #[test]
    fn test_bit_operations() {
        let mut m = BitMatrix::zeros(3, 3);

        // Test setting and getting bits
        m.set_bit(1, 2, true);
        assert!(m[(1, 2)]);
        assert!(!m[(1, 1)]);

        m.set_bit(0, 0, true);
        assert!(m[(0, 0)]);

        // Test setting bit to false
        m.set_bit(1, 2, false);
        assert!(!m[(1, 2)]);
    }

    #[test]
    fn test_is_zero() {
        let zero_matrix = BitMatrix::zeros(10, 10);
        assert!(zero_matrix.is_zero());

        let mut non_zero = BitMatrix::zeros(10, 10);
        non_zero.set_bit(5, 5, true);
        assert!(!non_zero.is_zero());

        // Test with identity matrix
        let identity = BitMatrix::identity(5);
        assert!(!identity.is_zero());
    }

    #[test]
    fn test_vstack_dimension_mismatch() {
        let m1 = BitMatrix::zeros(3, 4);
        let m2 = BitMatrix::zeros(2, 5); // Different number of columns
        m1.try_vstack(&m2).unwrap_err();
    }

    #[test]
    fn test_hstack_dimension_mismatch() {
        let m1 = BitMatrix::zeros(3, 4);
        let m2 = BitMatrix::zeros(5, 2); // Different number of rows
        m1.try_hstack(&m2).unwrap_err();
    }

    #[test]
    fn test_vstack_from_iter_empty() {
        let result = BitMatrix::vstack_from_iter(std::iter::empty());
        assert_eq!(result.rows(), 0);
        assert_eq!(result.cols(), 0);
    }

    #[test]
    fn test_hstack_from_iter_empty() {
        let result = BitMatrix::hstack_from_iter(std::iter::empty());
        assert_eq!(result.rows(), 0);
        assert_eq!(result.cols(), 0);
    }

    #[test]
    fn test_vstack_from_iter_single() {
        let m = BitMatrix::identity(3);
        let result = BitMatrix::vstack_from_iter([&m]);
        assert_eq!(result, m);
    }

    #[test]
    fn test_hstack_from_iter_single() {
        let m = BitMatrix::identity(3);
        let result = BitMatrix::hstack_from_iter([&m]);
        assert_eq!(result, m);
    }

    #[test]
    fn test_vstack_from_iter_multiple() {
        let m1 = BitMatrix::identity(2);
        let m2 = BitMatrix::zeros(3, 2);
        let m3 = BitMatrix::build(1, 2, |_, _| true);

        let result = BitMatrix::vstack_from_iter([&m1, &m2, &m3]);
        assert_eq!(result.rows(), 6);
        assert_eq!(result.cols(), 2);

        // Check that the stacking worked correctly
        assert_eq!(result[(0, 0)], true); // From identity
        assert_eq!(result[(1, 1)], true); // From identity
        assert_eq!(result[(2, 0)], false); // From zeros
        assert_eq!(result[(5, 1)], true); // From all-ones row
    }

    #[test]
    fn test_hstack_from_iter_multiple() {
        let m1 = BitMatrix::identity(2);
        let m2 = BitMatrix::zeros(2, 3);
        let m3 = BitMatrix::build(2, 1, |_, _| true);

        let result = BitMatrix::hstack_from_iter([&m1, &m2, &m3]);
        assert_eq!(result.rows(), 2);
        assert_eq!(result.cols(), 6);

        // Check that the stacking worked correctly
        assert_eq!(result[(0, 0)], true); // From identity
        assert_eq!(result[(1, 1)], true); // From identity
        assert_eq!(result[(0, 2)], false); // From zeros
        assert_eq!(result[(0, 5)], true); // From all-ones column
    }

    #[test]
    fn test_nullspace_empty_matrix() {
        let m = BitMatrix::zeros(0, 0);
        let nullspace = m.nullspace();
        assert_eq!(nullspace.len(), 0);
    }

    #[test]
    fn test_nullspace_zero_matrix() {
        let m = BitMatrix::zeros(3, 5);
        let nullspace = m.nullspace();
        assert_eq!(nullspace.len(), 5); // All columns are free variables

        // Verify each basis vector is in the nullspace
        for basis_vec in &nullspace {
            let result = &m * &basis_vec.transposed();
            assert!(result.is_zero());
        }
    }

    #[test]
    fn test_nullspace_identity_matrix() {
        let m = BitMatrix::identity(5);
        let nullspace = m.nullspace();
        assert_eq!(nullspace.len(), 0); // No nullspace for invertible matrix
    }

    #[test]
    fn test_nullspace_properties() {
        let mut rng = SmallRng::seed_from_u64(123);
        let m = BitMatrix::random(&mut rng, 4, 7);
        let nullspace = m.nullspace();

        // Each basis vector should be in the nullspace
        for basis_vec in &nullspace {
            let result = &m * &basis_vec.transposed();
            assert!(result.is_zero());
        }

        // The nullspace basis should be linearly independent
        if !nullspace.is_empty() {
            let ns_matrix = BitMatrix::vstack_from_iter(&nullspace);
            assert_eq!(ns_matrix.rank(), ns_matrix.rows());
        }
    }

    #[test]
    fn test_row_operations() {
        let mut m = BitMatrix::identity(3);
        let original = m.clone();

        // Test row swap
        m.swap_rows(0, 2);
        assert_eq!(m[(0, 0)], false);
        assert_eq!(m[(2, 2)], false);
        assert_eq!(m[(0, 2)], true);
        assert_eq!(m[(2, 0)], true);

        // Swap back
        m.swap_rows(0, 2);
        assert_eq!(m, original);

        // Test row addition (XOR)
        m.add_row(0, 1); // Add row 0 to row 1
        assert_eq!(m[(1, 0)], true); // XOR of 0 and 1
        assert_eq!(m[(1, 1)], true); // XOR of 1 and 0
    }

    #[test]
    fn test_add_bits_to_row() {
        let mut m = BitMatrix::zeros(3, 4);
        let bits = BitMatrix::build(1, 4, |_, j| j % 2 == 0);

        m.add_bits_to_row(bits.row(0), 1);

        // Check that the bits were added correctly
        for j in 0..4 {
            assert_eq!(m[(1, j)], j % 2 == 0);
        }
    }

    #[test]
    fn test_row_accessors() {
        let m = BitMatrix::identity(3);

        // Test immutable row access
        let row0 = m.row(0);
        assert_eq!(row0.len(), m.col_blocks);

        // Test mutable row access
        let mut m_mut = m.clone();
        {
            let row1 = m_mut.row_mut(1);
            // Modify the row (this is at the BitVec level)
            if row1.len() > 0 {
                row1[0] ^= MSB_ON; // Flip the first bit
            }
        }
        // The matrix should be modified
        assert_ne!(m_mut[(1, 0)], m[(1, 0)]);
    }

    #[test]
    fn test_transpose_inplace_rectangular_matrices() {
        let mut rng = SmallRng::seed_from_u64(665544);

        // Test various rectangular matrix dimensions
        let test_cases = [
            (1, 5),
            (5, 1),
            (3, 7),
            (7, 3),
            (10, 20),
            (20, 10),
            (32, 64),
            (64, 32),
            (63, 65),
            (65, 63),
            (100, 200),
            (200, 100),
            (128, 256),
            (256, 128),
        ];

        for (rows, cols) in test_cases {
            let original = BitMatrix::random(&mut rng, rows, cols);
            let expected = original.transposed();
            let mut actual = original.clone();
            actual.transpose_inplace();

            assert_eq!(actual.rows(), cols);
            assert_eq!(actual.cols(), rows);
            assert_eq!(actual, expected, "Failed for {}x{} matrix", rows, cols);

            // Double transpose should give original
            actual.transpose_inplace();
            assert_eq!(
                actual, original,
                "Double transpose failed for {}x{}",
                rows, cols
            );
        }
    }
}
