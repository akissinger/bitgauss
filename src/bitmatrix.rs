pub mod parallel;

use crate::bitvec::*;
use rand::Rng;
use std::{
    fmt,
    ops::{Index, Mul},
};

/// A matrix of bits, represented as a vector of blocks of bits
///
/// The matrix is stored in row-major order, with each row represented as a `BitRange` of `BitBlock`s. If
/// the number of columns is not a multiple of `BLOCKSIZE`, the last block in each row will be padded with 0s.
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

pub trait RowOps {
    fn add_row(&mut self, from: usize, to: usize);
    fn swap_rows(&mut self, from: usize, to: usize);
}

impl BitMatrix {
    #[inline]
    pub fn bit(&self, i: usize, j: usize) -> bool {
        self.data.bit(self.col_blocks * BLOCKSIZE * i + j)
    }

    #[inline]
    pub fn set_bit(&mut self, i: usize, j: usize, b: bool) {
        self.data.set_bit(self.col_blocks * BLOCKSIZE * i + j, b);
    }

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

    pub fn zeros(rows: usize, cols: usize) -> Self {
        let col_blocks = min_blocks(cols);
        BitMatrix {
            rows,
            cols,
            col_blocks,
            data: BitVec::zeros(rows * col_blocks),
        }
    }

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

    #[inline]
    pub fn random(rng: &mut impl Rng, rows: usize, cols: usize) -> Self {
        let col_blocks = min_blocks(cols);
        let num_blocks = BLOCKSIZE * rows * col_blocks;
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

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    #[inline]
    pub fn add_bits_to_row(&mut self, bits: &BitRange, row: usize) {
        self.data.xor_in(bits, row * self.col_blocks);
    }

    #[inline]
    pub fn row(&self, row: usize) -> &BitRange {
        &self.data[row * self.col_blocks..(row + 1) * self.col_blocks]
    }

    #[inline]
    pub fn row_mut(&mut self, row: usize) -> &mut BitRange {
        &mut self.data[row * self.col_blocks..(row + 1) * self.col_blocks]
    }

    #[inline]
    pub fn pad_to_square(&mut self) {
        let data_rows = self.data.len() / self.col_blocks;
        let row_blocks = min_blocks(data_rows);
        if data_rows != row_blocks * BLOCKSIZE || row_blocks != self.col_blocks {
            let blocks = usize::max(row_blocks, self.col_blocks);
            let mut data = Vec::with_capacity(BLOCKSIZE * blocks * blocks);
            for i in 0..(BLOCKSIZE * blocks) {
                for j in 0..blocks {
                    data.push(if i < self.rows() && j < self.col_blocks {
                        self.data[i * self.col_blocks + j]
                    } else {
                        0
                    });
                }
            }

            self.data = data.into();
            self.col_blocks = blocks;
        }
    }

    /// Main working function for transposition
    ///
    /// If `source` is given as `Some(bit_matrix)``, then copy bits from `bit_matrix` in transposed
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
                        buffer[k] = if i < m.data.len() { m.data[l] } else { 0 };
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

    /// Transpose the matrix in place, padding allocated memory with 0s if necessary
    #[inline]
    pub fn transpose_inplace(&mut self) {
        self.pad_to_square();
        self.transpose_helper(None);
    }

    /// Perform gaussian elimination while also performing matching row operations on `proxy`
    /// and returning a vector of pivot columns
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

    /// Perform gaussian elimination
    ///
    /// If `full` is true, then perform full Gauss-Jordan to produce reduced echelon form, otherwise
    /// just return echelon form
    #[inline]
    pub fn gauss(&mut self, full: bool) {
        self.gauss_helper(full, &mut ());
    }

    /// Compute the rank of the matrix using gaussian elimination
    #[inline]
    pub fn rank(&self) -> usize {
        self.clone().gauss_helper(false, &mut ()).len()
    }

    /// Compute the inverse of an invertible matrix
    pub fn inverse(&self) -> Self {
        if self.rows() != self.cols() {
            panic!("Matrix must be square");
        }
        let mut inv = BitMatrix::identity(self.cols());
        let pcols = self.clone().gauss_helper(true, &mut inv);

        if pcols.len() != self.cols() {
            panic!("Matrix is not invertible");
        }

        inv
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

        return true;
    }
}

impl Eq for BitMatrix {}

impl RowOps for () {
    #[inline]
    fn add_row(&mut self, _: usize, _: usize) {}

    #[inline]
    fn swap_rows(&mut self, _: usize, _: usize) {}
}

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

impl fmt::Display for BitMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, " {} ", if self[(i, j)] { 1 } else { 0 })?;
            }
            write!(f, "\n")?;
        }

        Ok(())
    }
}

impl Mul for &BitMatrix {
    type Output = BitMatrix;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.cols != rhs.rows {
            panic!(
                "Attempting to multiply matrices of incompatible dimensions: {} != {}",
                self.cols, rhs.rows
            );
        }
        let mut res = BitMatrix::zeros(self.rows, rhs.cols);

        for i in 0..self.rows {
            let row = res.row_mut(i);
            self.row(i).iter().enumerate().for_each(|(j, b)| {
                if b {
                    *row ^= rhs.row(j);
                }
            });
        }

        res
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{SeedableRng, rngs::SmallRng};

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
    fn pad_to_square() {
        let mut rng = SmallRng::seed_from_u64(1);
        let m = BitMatrix::random(&mut rng, 300, 200);
        let mut n = m.clone();
        n.pad_to_square();
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
}
