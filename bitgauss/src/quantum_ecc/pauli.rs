#![allow(clippy::suspicious_op_assign_impl)]
use std::ops::{AddAssign, MulAssign};

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum PauliLetter {
    I,
    X,
    Y,
    Z,
}

#[allow(dead_code)]
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum IPower {
    One,
    PlusI,
    NegOne,
    NegI,
}

impl IPower {
    fn to_u8(self) -> u8 {
        match self {
            IPower::One => 0,
            IPower::PlusI => 1,
            IPower::NegOne => 2,
            IPower::NegI => 3,
        }
    }
}

impl MulAssign<IPower> for IPower {
    fn mul_assign(&mut self, rhs: IPower) {
        *self += rhs.to_u8();
    }
}

impl AddAssign<u8> for IPower {
    #[inline]
    fn add_assign(&mut self, rhs: u8) {
        match rhs % 4 {
            0 => {}
            1 => match self {
                IPower::One => {
                    *self = IPower::PlusI;
                }
                IPower::PlusI => {
                    *self = IPower::NegOne;
                }
                IPower::NegOne => {
                    *self = IPower::NegI;
                }
                IPower::NegI => {
                    *self = IPower::One;
                }
            },
            2 => match self {
                IPower::One => {
                    *self = IPower::NegOne;
                }
                IPower::PlusI => {
                    *self = IPower::NegI;
                }
                IPower::NegOne => {
                    *self = IPower::One;
                }
                IPower::NegI => {
                    *self = IPower::PlusI;
                }
            },
            3 => match self {
                IPower::One => {
                    *self = IPower::NegI;
                }
                IPower::PlusI => {
                    *self = IPower::One;
                }
                IPower::NegOne => {
                    *self = IPower::PlusI;
                }
                IPower::NegI => {
                    *self = IPower::NegOne;
                }
            },
            _ => unreachable!(),
        }
    }
}
