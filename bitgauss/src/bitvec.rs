use rand::Rng;
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
/// use bitgauss::BitVec;
///
/// // Create a BitVec of 256 bits, all set to zero
/// let mut bv = BitVec::zeros(4);
/// bv.set_bit(5, true);
/// assert!(bv.bit(5));
/// ```
///
/// # Note
///
/// Many methods are implemented via dereferencing to [`BitSlice`], which provides
/// additional bitwise and range operations.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct BitVec(Vec<BitBlock>);

/// A range of bits, represented as a slice of [`BitBlock`]s.
///
/// Provides methods for bitwise operations, iteration, and bit access within the range.
#[derive(RefCast, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(transparent)]
pub struct BitSlice([BitBlock]);

/// Iterator over the bits in a [`BitSlice`].
///
/// Yields each bit as a `bool`, starting from the most significant bit of the first block.
pub struct BitIter<'a> {
    inner: std::slice::Iter<'a, BitBlock>,
    c: usize,
    block: BitBlock,
}
impl<'a> Iterator for BitIter<'a> {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        if self.c == BLOCKSIZE {
            self.block = self.inner.next().copied()?;
            self.c = 0;
        }
        let bit = self.block & MSB_ON == MSB_ON;
        self.block = self.block.wrapping_shl(1);
        self.c += 1;
        Some(bit)
    }
}

pub type BitBlockIter<'a> = std::iter::Copied<std::slice::Iter<'a, BitBlock>>;

impl BitSlice {
    /// Returns a copy of the range as a [`BitVec`].
    #[inline]
    pub fn to_vec(&self) -> BitVec {
        self.0.to_vec().into()
    }

    /// Returns an iterator over the [`BitBlock`]s in this range.
    #[inline]
    pub fn block_iter(&self) -> BitBlockIter<'_> {
        self.0.iter().copied()
    }

    /// Returns a mutable iterator over the [`BitBlock`]s in this range.
    #[inline]
    pub fn block_iter_mut(&mut self) -> impl Iterator<Item = &mut BitBlock> {
        self.0.iter_mut()
    }

    /// Returns an iterator over all bits in this range as `bool`s.
    #[inline]
    pub fn iter(&self) -> BitIter {
        BitIter {
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

    /// Computes the dot product (mod 2) of two [`BitSlice`]s.
    ///
    /// Returns `true` if the number of matching 1s is odd, otherwise `false`.
    #[inline]
    pub fn dot(&self, rhs: &BitSlice) -> bool {
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

    /// Performs an XOR operation between source and target ranges.
    pub fn xor_range(&mut self, source: usize, target: usize, len: usize) {
        for i in 0..len {
            self.0[target + i] ^= self.0[source + i];
        }
    }

    /// Extracts a subrange of bit blocks into a new [`BitVec`].
    pub fn extract(&self, start: usize, len: usize) -> BitVec {
        BitVec(self.0[start..(start + len)].into())
    }

    /// XORs another [`BitSlice`] into self, starting at a given target position.
    pub fn xor_in(&mut self, source: &BitSlice, target_pos: usize) {
        for i in 0..source.len() {
            self.0[target_pos + i] ^= source.0[i];
        }
    }

    /// Swaps two bit blocks at given indices.
    #[inline]
    pub fn swap(&mut self, source: usize, target: usize) {
        self.0.swap(source, target);
    }

    /// Swaps ranges of bit blocks.
    #[inline]
    pub fn swap_range(&mut self, source: usize, target: usize, len: usize) {
        for i in 0..len {
            self.0.swap(source + i, target + i);
        }
    }

    /// Returns the number of [`BitBlock`]s in the range.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns true if the range contains no blocks.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the total number of bits in the range.
    #[inline]
    pub fn num_bits(&self) -> usize {
        self.0.len() * BLOCKSIZE
    }
}

impl Index<Range<usize>> for BitSlice {
    type Output = BitSlice;
    fn index(&self, index: Range<usize>) -> &Self::Output {
        BitSlice::ref_cast(&self.0[index])
    }
}

impl Index<usize> for BitSlice {
    type Output = BitBlock;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}

impl IndexMut<usize> for BitSlice {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.0.index_mut(index)
    }
}

impl IndexMut<Range<usize>> for BitSlice {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        BitSlice::ref_cast_mut(self.0.index_mut(index))
    }
}

impl BitVec {
    /// Returns the number of [`BitBlock`]s in the vector.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the total number of bits in the vector.
    #[inline]
    pub fn num_bits(&self) -> usize {
        self.0.len() * BLOCKSIZE
    }

    /// Returns true if the vector contains no blocks.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value of the bit at the specified index.
    #[inline]
    pub fn bit(&self, index: usize) -> bool {
        self.deref().bit(index)
    }

    /// Sets the bit at the given index to the provided value.
    #[inline]
    pub fn set_bit(&mut self, index: usize, value: bool) {
        self.deref_mut().set_bit(index, value)
    }

    /// Returns an iterator over all bits in this vector as `bool`s.
    #[inline]
    pub fn iter(&self) -> BitIter {
        self.deref().iter()
    }

    /// Counts the number of bits set to 1 in the entire vector.
    #[inline]
    pub fn count_ones(&self) -> u32 {
        self.deref().count_ones()
    }

    /// Counts the number of bits set to 0 in the entire vector.
    #[inline]
    pub fn count_zeros(&self) -> u32 {
        self.deref().count_zeros()
    }

    /// Returns a [`BitSlice`] that represents a subrange of the vector.
    ///
    /// # Arguments
    /// * `from_block` - Starting block index (inclusive).
    /// * `to_block` - Ending block index (exclusive).
    #[inline]
    pub fn bit_range(&self, from_block: usize, to_block: usize) -> &BitSlice {
        BitSlice::ref_cast(&self.0[from_block..to_block])
    }

    /// Returns a mutable [`BitSlice`] that represents a subrange of the vector.
    ///
    /// # Arguments
    /// * `from_block` - Starting block index (inclusive).
    /// * `to_block` - Ending block index (exclusive).
    #[inline]
    pub fn bit_range_mut(&mut self, from_block: usize, to_block: usize) -> &mut BitSlice {
        BitSlice::ref_cast_mut(&mut self.0[from_block..to_block])
    }

    /// Constructs a random [`BitVec`] with the specified number of [`BitBlock`]s.
    ///
    /// # Arguments
    /// * `rng` - A mutable reference to a random number generator.
    /// * `num_blocks` - The block size of the generated bit vector.
    #[inline]
    pub fn random(rng: &mut impl Rng, num_blocks: usize) -> Self {
        (0..num_blocks).map(|_| rng.random::<BitBlock>()).collect()
    }

    /// Constructs a [`BitVec`] with all bits set to zero.
    ///
    /// # Arguments
    /// * `num_blocks` - The block size of the new bit vector.
    #[inline]
    pub fn zeros(num_blocks: usize) -> Self {
        BitVec(vec![0; num_blocks])
    }

    /// Constructs a [`BitVec`] with all bits set to one.
    ///
    /// # Arguments
    /// * `num_blocks` - The block size of the new bit vector.
    #[inline]
    pub fn ones(num_blocks: usize) -> Self {
        BitVec(vec![BitBlock::MAX; num_blocks])
    }

    /// Constructs a new empty [`BitVec`].
    pub fn new() -> Self {
        BitVec(Vec::new())
    }

    /// Constructs a new [`BitVec`] with the specified capacity in blocks
    pub fn with_capacity(num_blocks: usize) -> Self {
        BitVec(Vec::with_capacity(num_blocks))
    }

    /// Reserves capacity for at least `additional` more blocks in the vector
    pub fn reserve(&mut self, additional: usize) {
        self.0.reserve(additional);
    }

    /// Push a new block to the end of the vector
    pub fn push_block(&mut self, block: BitBlock) {
        self.0.push(block);
    }

    /// Extends a [`BitVec`] with the contents of a [`BitSlice`]
    pub fn extend_from_slice(&mut self, other: &BitSlice) {
        self.0.extend_from_slice(&other.0);
    }

    /// Extends a [`BitVec`] with the contents of a [`BitSlice`], left-shifting the bits in each block
    ///
    /// Note this method assumes that the last `shift` bits in `self` are zero
    pub fn extend_from_slice_left_shifted(&mut self, other: &BitSlice, shift: usize) {
        if shift >= BLOCKSIZE {
            panic!("Shift must be less than BLOCKSIZE");
        } else if shift == 0 {
            self.extend_from_slice(other);
            return;
        } else if self.0.is_empty() {
            panic!("Cannot append to an empty BitVec with left shift");
        }

        self.0.reserve(other.0.len());
        for bits in other.0.iter() {
            let left_part = bits.wrapping_shr((BLOCKSIZE - shift) as u32);
            let right_part = bits.wrapping_shl(shift as u32);
            if let Some(last) = self.0.last_mut() {
                *last |= left_part;
            }
            self.0.push(right_part);
        }
    }

    /// Pops the last block from the vector and returns it
    pub fn pop(&mut self) -> Option<BitBlock> {
        self.0.pop()
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

impl BitAndAssign<&Self> for BitSlice {
    #[inline]
    fn bitand_assign(&mut self, rhs: &Self) {
        for (bits0, bits1) in self.0.iter_mut().zip(rhs.0.iter()) {
            *bits0 &= bits1;
        }
    }
}

impl BitXorAssign<&Self> for BitSlice {
    #[inline]
    fn bitxor_assign(&mut self, rhs: &BitSlice) {
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
    type Target = BitSlice;
    fn deref(&self) -> &Self::Target {
        BitSlice::ref_cast(&self.0)
    }
}

impl DerefMut for BitVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        BitSlice::ref_cast_mut(&mut self.0)
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
    use rand::{rngs::SmallRng, SeedableRng};

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
        // let r: &BitSlice = &vec;
        let r1: &BitSlice = &vec[4..9];

        for i in 0..r1.len() {
            assert_eq!(vec[4 + i], r1[i]);
        }
    }

    // test extend_from_slice_left_shifted
    #[test]
    fn extend_from_slice_left_shifted() {
        let mut rng = SmallRng::seed_from_u64(1);
        let shift = 17;
        let mask = BitBlock::MAX.wrapping_shl(17);

        let mut v1 = BitVec::random(&mut rng, 10);
        v1[9] &= mask;

        let v2 = BitVec::random(&mut rng, 10);

        let mut v3 = v1.clone();
        v3.extend_from_slice_left_shifted(&v2, shift);

        for i in 0..v3.num_bits() {
            if i < 10 * BLOCKSIZE - shift {
                assert_eq!(v3.bit(i), v1.bit(i));
            } else if i < 20 * BLOCKSIZE - shift {
                assert_eq!(v3.bit(i), v2.bit(i - (10 * BLOCKSIZE - shift)));
            } else {
                assert_eq!(v3.bit(i), false);
            }
        }
    }
}
