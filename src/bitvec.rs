use rand::Rng;
use rayon::prelude::*;
use ref_cast::RefCast;
use std::fmt;
pub use std::ops::{BitAndAssign, BitXorAssign, Deref, DerefMut, Index, IndexMut, Range};

/// A block of bits. This is an alias for [`u64`]
pub type BitBlock = u64;

/// Number of bits in a [`BitBlock`]
pub const BLOCKSIZE: usize = 64;

/// Bitwise AND with this constant to set most signficant bit to zero
pub const MSB_OFF: BitBlock = 0x7fffffffffffffff;

/// Bitwise OR with this constant to set most signficant bit to one
pub const MSB_ON: BitBlock = 0x8000000000000000;

/// Returns the minimum number of [`BitBlock`]s required to store the given number of bits.
///
/// # Arguments
///
/// * `bits` - The number of bits to store.
///
/// # Returns
///
/// The minimum number of [`BitBlock`]s (each of size [`BLOCKSIZE`]) needed to store `bits` bits.
/// If `bits` is not a multiple of [`BLOCKSIZE`], the result is rounded up to ensure all bits fit.
#[inline]
pub fn min_blocks(bits: usize) -> usize {
    bits / BLOCKSIZE + if bits % BLOCKSIZE == 0 { 0 } else { 1 }
}

/// A vector of bits, stored efficiently as a vector of [`BitBlock`]s (which alias to `u64`).
///
/// `BitVec` provides a compact and performant way to store and manipulate large bit vectors.
/// It supports bitwise operations, random and zero/one initialization, and conversion to and from
/// boolean vectors. The bits are packed into 64-bit blocks, and the struct offers methods for
/// accessing, setting, and iterating over individual bits or ranges of bits.
///
/// # Examples
///
/// ```
/// use bitgauss::bitvec::*;
///
/// // Create a BitVec of 256 bits, all set to zero
/// let mut bv = BitVec::zeros(4);
/// bv.set_bit(5, true);
/// assert!(bv.bit(5));
/// ```
///
/// # Note
///
/// Many methods are implemented via dereferencing to [`BitRange`], which provides
/// additional bitwise and range operations.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct BitVec(Vec<BitBlock>);

/// A range of bits, represented as a slice of [`BitBlock`]s.
///
/// Provides methods for bitwise operations, iteration, and bit access within the range.
#[derive(RefCast, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(transparent)]
pub struct BitRange([BitBlock]);

/// Iterator over the bits in a [`BitRange`].
///
/// Yields each bit as a `bool`, starting from the most significant bit of the first block.
pub struct BitRangeIter<'a> {
    inner: std::slice::Iter<'a, BitBlock>,
    c: usize,
    block: BitBlock,
}
impl<'a> Iterator for BitRangeIter<'a> {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        if self.c == BLOCKSIZE {
            self.block = self.inner.next().copied()?;
            self.c = 0;
        }
        let bit = self.block & MSB_ON == MSB_ON;
        self.block <<= 1;
        self.c += 1;
        Some(bit)
    }
}

impl BitRange {
    /// Returns a copy of the range as a [`BitVec`].
    #[inline]
    pub fn to_vec(&self) -> BitVec {
        self.0.to_vec().into()
    }

    /// Divides the range into mutable parallel chunks of the given size.
    ///
    /// Useful for parallel processing over disjoint bit regions.
    ///
    /// # Arguments
    ///
    /// * `chunk_size` - Number of blocks per chunk.
    #[inline]
    pub fn par_chunks_mut(
        &mut self,
        chunk_size: usize,
    ) -> impl ParallelIterator<Item = &mut BitRange> {
        self.0
            .par_chunks_mut(chunk_size)
            .map(|x| BitRange::ref_cast_mut(x))
    }

    /// Returns an iterator over the [`BitBlock`]s in this range.
    #[inline]
    pub fn block_iter(&self) -> impl Iterator<Item = BitBlock> {
        self.0.iter().copied()
    }

    /// Returns a mutable iterator over the [`BitBlock`]s in this range.
    #[inline]
    pub fn block_iter_mut(&mut self) -> impl Iterator<Item = &mut BitBlock> {
        self.0.iter_mut()
    }

    /// Returns an iterator over all bits in this range as `bool`s.
    #[inline]
    pub fn iter(&self) -> BitRangeIter {
        BitRangeIter {
            inner: self.0.iter(),
            c: BLOCKSIZE,
            block: 0,
        }
    }

    /// Counts the number of bits set to 1 in the entire range.
    #[inline]
    pub fn count_ones(&self) -> u32 {
        self.block_iter().fold(0, |c, bits| c + bits.count_ones())
    }

    /// Counts the number of bits set to 0 in the entire range.
    #[inline]
    pub fn count_zeros(&self) -> u32 {
        self.block_iter().fold(0, |c, bits| c + bits.count_zeros())
    }

    /// Computes the dot product (mod 2) of two [`BitRange`]s.
    ///
    /// Returns `true` if the number of matching 1s is odd, otherwise `false`.
    #[inline]
    pub fn dot(&self, rhs: &BitRange) -> bool {
        let mut c = 0;
        for (bits0, bits1) in self.0.iter().zip(rhs.0.iter()) {
            c ^= (*bits0 & *bits1).count_ones() & 1;
        }

        c == 1
    }

    /// Returns the value of the bit at the specified index.
    ///
    /// # Panics
    ///
    /// Panics if the index is out of range.
    #[inline]
    pub fn bit(&self, index: usize) -> bool {
        let block_index = index / BLOCKSIZE;
        let bit_index = (index % BLOCKSIZE) as u32;
        let block = self.0[block_index].rotate_left(bit_index);
        block & MSB_ON == MSB_ON
    }

    /// Sets the bit at the given index to the provided value.
    ///
    /// # Arguments
    ///
    /// * `index` - Bit index to set.
    /// * `value` - `true` to set to 1, `false` to set to 0.
    ///
    /// # Panics
    ///
    /// Panics if the index is out of range.
    #[inline]
    pub fn set_bit(&mut self, index: usize, value: bool) {
        let block_index = index / BLOCKSIZE;
        let bit_index = (index % BLOCKSIZE) as u32;
        let mut block = self.0[block_index].rotate_left(bit_index);
        if value {
            block |= MSB_ON;
        } else {
            block &= MSB_OFF;
        }
        self.0[block_index] = block.rotate_right(bit_index);
    }

    /// Returns the position (in bits) of the first 1-bit in the specified range of [`BitBlock`]s.
    ///
    /// # Arguments
    ///
    /// * `from` - Starting block index.
    /// * `to` - Ending block index (exclusive).
    ///
    /// # Returns
    ///
    /// `Some(bit_index)` if a 1-bit is found, otherwise `None`.
    pub fn first_one_in_range(&self, from: usize, to: usize) -> Option<usize> {
        for i in from..to {
            if self.0[i] != 0 {
                return Some((i - from) * BLOCKSIZE + (self.0[i].leading_zeros() as usize));
            }
        }
        None
    }

    pub fn xor_range(&mut self, source: usize, target: usize, len: usize) {
        for i in 0..len {
            self.0[target + i] ^= self.0[source + i];
        }
    }

    pub fn extract(&self, start: usize, len: usize) -> BitVec {
        BitVec(self.0[start..(start + len)].into())
    }

    pub fn xor_in(&mut self, source: &BitRange, target_pos: usize) {
        for i in 0..source.len() {
            self.0[target_pos + i] ^= source.0[i];
        }
    }

    #[inline]
    pub fn swap(&mut self, source: usize, target: usize) {
        self.0.swap(source, target);
    }

    #[inline]
    pub fn swap_range(&mut self, source: usize, target: usize, len: usize) {
        for i in 0..len {
            self.0.swap(source + i, target + i);
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn num_bits(&self) -> usize {
        self.0.len() * BLOCKSIZE
    }
}

impl Index<Range<usize>> for BitRange {
    type Output = BitRange;
    fn index(&self, index: Range<usize>) -> &Self::Output {
        BitRange::ref_cast(&self.0[index])
    }
}

impl Index<usize> for BitRange {
    type Output = BitBlock;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}

impl IndexMut<usize> for BitRange {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.0.index_mut(index)
    }
}

impl IndexMut<Range<usize>> for BitRange {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        BitRange::ref_cast_mut(self.0.index_mut(index))
    }
}

impl BitVec {
    #[inline]
    pub fn bit_range(&self, from_block: usize, to_block: usize) -> &BitRange {
        BitRange::ref_cast(&self.0[from_block..to_block])
    }

    #[inline]
    pub fn bit_range_mut(&mut self, from_block: usize, to_block: usize) -> &mut BitRange {
        BitRange::ref_cast_mut(&mut self.0[from_block..to_block])
    }

    #[inline]
    pub fn random(rng: &mut impl Rng, num_blocks: usize) -> Self {
        (0..num_blocks).map(|_| rng.random::<BitBlock>()).collect()
    }

    #[inline]
    pub fn zeros(num_blocks: usize) -> Self {
        BitVec(vec![0; num_blocks])
    }

    #[inline]
    pub fn ones(num_blocks: usize) -> Self {
        BitVec(vec![BitBlock::MAX; num_blocks])
    }
}

impl fmt::Display for BitVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &bits in self.0.iter() {
            write!(f, "{:064b}", bits)?;
        }
        Ok(())
    }
}

impl BitAndAssign<&Self> for BitRange {
    #[inline]
    fn bitand_assign(&mut self, rhs: &Self) {
        for (bits0, bits1) in self.0.iter_mut().zip(rhs.0.iter()) {
            *bits0 &= bits1;
        }
    }
}

impl BitXorAssign<&Self> for BitRange {
    #[inline]
    fn bitxor_assign(&mut self, rhs: &BitRange) {
        for (bits0, bits1) in self.0.iter_mut().zip(rhs.0.iter()) {
            *bits0 ^= bits1;
        }
    }
}

impl From<Vec<BitBlock>> for BitVec {
    fn from(value: Vec<BitBlock>) -> Self {
        BitVec(value)
    }
}

impl From<BitVec> for Vec<BitBlock> {
    fn from(value: BitVec) -> Self {
        value.0
    }
}

impl FromIterator<BitBlock> for BitVec {
    fn from_iter<T: IntoIterator<Item = BitBlock>>(iter: T) -> Self {
        Vec::from_iter(iter).into()
    }
}

impl FromIterator<bool> for BitVec {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let mut v = vec![];
        let mut c = 0;
        let mut block: BitBlock = 0;
        for bit in iter {
            if bit {
                block |= 1;
            }
            c += 1;
            if c == BLOCKSIZE {
                c = 0;
                v.push(block);
                block = 0;
            } else {
                block <<= 1;
            }
        }

        if c != 0 {
            block <<= BLOCKSIZE - c - 1;
            v.push(block);
        }

        BitVec(v)
    }
}

impl Deref for BitVec {
    type Target = BitRange;
    fn deref(&self) -> &Self::Target {
        BitRange::ref_cast(&self.0)
    }
}

impl DerefMut for BitVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        BitRange::ref_cast_mut(&mut self.0)
    }
}

impl From<Vec<bool>> for BitVec {
    fn from(value: Vec<bool>) -> Self {
        BitVec::from_iter(value.iter().copied())
    }
}

impl From<BitVec> for Vec<bool> {
    fn from(value: BitVec) -> Self {
        value.iter().collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{SeedableRng, rngs::SmallRng};

    #[test]
    fn bit_xor_and() {
        let sz = 8;
        let mut rng = SmallRng::seed_from_u64(1);
        let vec = BitVec::random(&mut rng, sz);
        let mut vec1 = vec.clone();
        *vec1 ^= &vec;
        assert_eq!(vec1, BitVec::zeros(sz));

        vec1 = vec.clone();
        *vec1 &= &BitVec::zeros(sz);
        assert_eq!(vec1, BitVec::zeros(sz));

        vec1 = vec.clone();
        *vec1 &= &vec;
        assert_eq!(vec1, vec);
    }

    #[test]
    fn bit_get_set() {
        let sz = 4;
        let bits = vec![0, 3, 100, 201, 255];

        let mut vec0 = BitVec::zeros(sz);
        for &b in &bits {
            vec0.set_bit(b, true);
        }

        for i in 0..(sz * BLOCKSIZE) {
            assert_eq!(vec0.bit(i), bits.contains(&i));
        }

        let mut vec1 = BitVec::ones(sz);
        for &b in &bits {
            vec1.set_bit(b, false);
        }

        for i in 0..(sz * BLOCKSIZE) {
            assert_eq!(vec1.bit(i), !bits.contains(&i));
        }
    }

    #[test]
    fn bool_vec() {
        let mut rng = SmallRng::seed_from_u64(1);
        let bool_vec: Vec<bool> = (0..300).map(|_| rng.random()).collect();
        let vec: BitVec = bool_vec.clone().into();
        let bool_vec1: Vec<bool> = vec.clone().into();

        // converting to BitVec will pad to a multiple of BLOCKSIZE
        for (i, &b) in bool_vec.iter().enumerate() {
            assert_eq!((i, vec.bit(i)), (i, b));
            assert_eq!((i, bool_vec1[i]), (i, b));
        }

        // ...so the remaining bits should be 0
        assert_eq!(vec.num_bits(), bool_vec1.len());
        for i in bool_vec.len()..vec.len() {
            assert_eq!((i, vec.bit(i)), (i, false));
            assert_eq!((i, bool_vec1[i]), (i, false));
        }
    }

    #[test]
    fn xor_range() {
        let i = BitBlock::MAX;
        let vec0: BitVec = vec![0, i, 0, i, 0, 0, i, i, 0, 0].into();

        let mut vec1 = vec0.clone();
        vec1.xor_range(1, 5, 3);

        let vec2: BitVec = vec![0, i, 0, i, 0, i, i, 0, 0, 0].into();
        assert_eq!(vec1, vec2);

        vec1.xor_range(1, 5, 3);
        assert_eq!(vec0, vec1);
    }

    #[test]
    fn block_index() {
        let mut rng = SmallRng::seed_from_u64(1);
        let vec: BitVec = BitVec::random(&mut rng, 10);
        // let r: &BitRange = &vec;
        let r1: &BitRange = &vec[4..9];

        for i in 0..r1.len() {
            assert_eq!(vec[4 + i], r1[i]);
        }
    }
}
