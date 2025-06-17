//! This crate provides two main data structures:
//! - [`crate::bitvec::BitVec`]: a vector of bits stored in 64-bit chunks, along with convenience
//!   methods for indexing, slicing, and manipulating bits
//! - [`crate::bitmatrix::BitMatrix`]: a two-dimensional matrix based on `BitVec`, which implements
//!   basic linear algebraic operations

pub mod bitmatrix;
pub mod bitvec;
