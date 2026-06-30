use std::{
    fmt::Display,
    ops::{Add, AddAssign},
};

use crate::data::{min_blocks, BitData};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct F2nSymplectic {
    pub(crate) first_part: BitData,
    pub(crate) second_part: BitData,
    pub(crate) my_n: usize,
}

impl Display for F2nSymplectic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "Xs: {:?}",
            self.first_part
                .into_iter()
                .take(self.my_n)
                .map(|b| if b { '1' } else { '0' })
                .collect::<String>()
        ))?;
        f.write_fmt(format_args!(
            "Zs: {:?}",
            self.second_part
                .into_iter()
                .take(self.my_n)
                .map(|b| if b { '1' } else { '0' })
                .collect::<String>()
        ))?;
        Ok(())
    }
}

impl F2nSymplectic {
    /// The standard form on `F_2^{2n}`
    ///
    /// # Panics
    ///
    /// The two must be in the same vector space with the same `n`.
    pub fn omega(&self, other: &Self) -> bool {
        assert_eq!(self.my_n, other.my_n);
        let a = self.first_part.dot(&other.second_part);
        let b = self.second_part.dot(&other.first_part);
        a ^ b
    }

    pub(crate) fn weight(&self) -> usize {
        let mut either_ones = self.first_part.clone();
        *either_ones |= &self.second_part;
        either_ones.count_ones()
    }

    /// Conjugate by `H` on the specified qubits.
    /// If on all of them, we do them all at once with a single `mem::swap`.
    /// If on none of them, we do nothing.
    /// Otherwise we do the switching of roles of `X` and `Z` only for that qubit.
    pub fn hadamard_transform(&mut self, mut where_hs: BitData) {
        if where_hs.len() < self.my_n {
            where_hs.extend_from_slice(&BitData::zeros(min_blocks(self.my_n)));
        }
        if where_hs.count_zeros() == 0 {
            core::mem::swap(&mut self.first_part, &mut self.second_part);
            return;
        }
        if where_hs.count_ones() == 0 {
            return;
        }
        for (idx, cur_h_bit) in where_hs.into_iter().take(self.my_n).enumerate() {
            if cur_h_bit {
                let from_first = self.first_part.bit(idx);
                let from_second = self.second_part.bit(idx);
                self.first_part.set_bit(idx, from_second);
                self.second_part.set_bit(idx, from_first);
            }
        }
    }

    /// The inclusion of `F_2^{2n} \to F_2^{2(n+r)}` with `0`s in the last `r`
    /// parts for both the `X` and `Z` parts which are the two standard Lagrangian
    /// subspaces.
    #[allow(clippy::missing_panics_doc)]
    pub fn pad_right(&mut self, mut right_pad: usize) {
        assert!(
            self.first_part.num_bits() >= self.my_n,
            "{} vs {}",
            self.first_part.num_bits(),
            self.my_n
        );
        let how_much_room_already = self.first_part.num_bits().saturating_sub(self.my_n);
        if right_pad > how_much_room_already {
            self.my_n += right_pad;
            right_pad = right_pad.saturating_sub(how_much_room_already);
            let requisite_zeros = BitData::zeros(min_blocks(right_pad));
            self.first_part.extend_from_slice(&requisite_zeros);
            self.second_part.extend_from_slice(&requisite_zeros);
        } else {
            self.my_n += right_pad;
        }
    }
}

impl Add for &F2nSymplectic {
    type Output = F2nSymplectic;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.my_n, rhs.my_n);
        let mut first_part = self.first_part.clone();
        let mut second_part = self.second_part.clone();
        *first_part ^= &rhs.first_part;
        *second_part ^= &rhs.second_part;
        F2nSymplectic {
            first_part,
            second_part,
            my_n: self.my_n,
        }
    }
}

impl AddAssign<&Self> for F2nSymplectic {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.my_n, rhs.my_n);
        *self.first_part ^= &rhs.first_part;
        *self.second_part ^= &rhs.second_part;
    }
}
