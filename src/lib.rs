//! `bitgauss` is a Rust library for doing linear algebra over the 2-element finite field. Some features include:
//! - getting and setting individual matrix elements (as `bool`s)
//! - fast row operations and dot product using bitwise operations
//! - fast in-place and out-of-place matrix transpose using a [recursive block method](https://github.com/dsnet/matrix-transpose)
//! - horizontal and vertical concatenation of matrices
//! - matrix multiplication
//! - Gaussian elimination and related methods (e.g. rank and inverse)
//!
//! The two main data structures provided by this crate are:
//! - [`BitVec`]: a vector of bits stored in 64-bit chunks, along with convenience
//!   methods for indexing, slicing, and manipulating bits
//! - [`BitMatrix`]: a two-dimensional matrix based on `BitVec`, which implements
//!   basic linear algebraic operations

#[allow(clippy::needless_range_loop, clippy::suspicious_arithmetic_impl)]
pub mod bitmatrix;
pub mod bitvec;

pub use bitmatrix::{BitMatrix, RowOps};
pub use bitvec::{BitBlock, BitSlice, BitVec};
