use std::{fmt::Display, ops::MulAssign};

use crate::{
    data::{min_blocks, BitData},
    quantum_ecc::{
        f2n_symplectic::F2nSymplectic,
        pauli::{IPower, PauliLetter},
    },
    BitMatrix, RowOps,
};

impl From<&F2nSymplectic> for Vec<PauliLetter> {
    fn from(f2n_symplectic: &F2nSymplectic) -> Self {
        let mut to_return = Vec::with_capacity(f2n_symplectic.my_n);
        let first_iter = f2n_symplectic
            .first_part
            .into_iter()
            .take(f2n_symplectic.my_n);
        let second_iter = f2n_symplectic
            .second_part
            .into_iter()
            .take(f2n_symplectic.my_n);
        for (x_part, z_part) in first_iter.zip(second_iter) {
            match (x_part, z_part) {
                (true, true) => {
                    to_return.push(PauliLetter::Y);
                }
                (true, false) => {
                    to_return.push(PauliLetter::X);
                }
                (false, true) => {
                    to_return.push(PauliLetter::Z);
                }
                (false, false) => {
                    to_return.push(PauliLetter::I);
                }
            }
        }
        to_return
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct PauliString {
    #[allow(dead_code)]
    overall_factor: IPower,
    pauli_letters: F2nSymplectic,
}

impl Display for PauliString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{:?} * ", self.overall_factor))?;
        f.write_fmt(format_args!(
            "{:?}",
            Vec::<PauliLetter>::from(&self.pauli_letters)
        ))?;
        Ok(())
    }
}

impl PauliString {
    /// Parse a shorthand string like XZZIXX
    /// with no prefactor.
    /// If you want the prefactor, parse and then multiply by `IPower` after that
    pub fn parse(s: &str) -> Option<Self> {
        let mut pauli_letters = Vec::with_capacity(s.len());
        for letter in s.chars() {
            if letter.is_whitespace() || letter.is_ascii_punctuation() {
            } else if letter == '0' || letter == 'I' || letter == 'i' {
                pauli_letters.push(PauliLetter::I);
            } else if letter == 'X' || letter == 'x' {
                pauli_letters.push(PauliLetter::X);
            } else if letter == 'Y' || letter == 'y' {
                pauli_letters.push(PauliLetter::Y);
            } else if letter == 'Z' || letter == 'z' {
                pauli_letters.push(PauliLetter::Z);
            } else {
                return None;
            }
        }
        Some(Self {
            overall_factor: IPower::One,
            pauli_letters: F2nSymplectic::from(pauli_letters.into_iter()),
        })
    }

    pub fn weight_one(
        where_nontrivial: usize,
        how_many_total: usize,
        which_letter: PauliLetter,
        overall_factor: IPower,
    ) -> Self {
        let mut x_vector = BitData::zeros(min_blocks(how_many_total));
        let mut z_vector = BitData::zeros(min_blocks(how_many_total));
        match which_letter {
            PauliLetter::I => {}
            PauliLetter::X => {
                x_vector.set_bit(where_nontrivial, true);
            }
            PauliLetter::Y => {
                x_vector.set_bit(where_nontrivial, true);
                z_vector.set_bit(where_nontrivial, true);
            }
            PauliLetter::Z => {
                z_vector.set_bit(where_nontrivial, true);
            }
        }
        Self {
            overall_factor,
            pauli_letters: F2nSymplectic {
                first_part: x_vector,
                second_part: z_vector,
                my_n: how_many_total,
            },
        }
    }

    /// Do these two commute in `P_n`
    pub fn commutes_with(&self, other: &Self) -> bool {
        !self.pauli_letters.omega(&other.pauli_letters)
    }

    pub fn weight(&self) -> usize {
        self.pauli_letters.weight()
    }
}

impl<T> From<T> for F2nSymplectic
where
    T: ExactSizeIterator<Item = PauliLetter>,
{
    fn from(pauli_letters: T) -> Self {
        let num_qubits = pauli_letters.len();
        let mut x_part = Vec::with_capacity(num_qubits);
        let mut z_part = Vec::with_capacity(num_qubits);
        for cur_letter in pauli_letters {
            match cur_letter {
                PauliLetter::I => {
                    x_part.push(false);
                    z_part.push(false);
                }
                PauliLetter::X => {
                    x_part.push(true);
                    z_part.push(false);
                }
                PauliLetter::Y => {
                    x_part.push(true);
                    z_part.push(true);
                }
                PauliLetter::Z => {
                    x_part.push(false);
                    z_part.push(true);
                }
            }
        }
        Self {
            first_part: x_part.into(),
            second_part: z_part.into(),
            my_n: num_qubits,
        }
    }
}

impl MulAssign<IPower> for PauliString {
    fn mul_assign(&mut self, rhs: IPower) {
        self.overall_factor *= rhs;
    }
}

impl MulAssign<&PauliString> for PauliString {
    fn mul_assign(&mut self, rhs: &PauliString) {
        let mut self_is_y = self.pauli_letters.first_part.clone();
        *self_is_y &= &self.pauli_letters.second_part;
        let mut self_is_x = self.pauli_letters.second_part.negate();
        *self_is_x &= &self.pauli_letters.first_part;
        let mut self_is_z = self.pauli_letters.first_part.negate();
        *self_is_z &= &self.pauli_letters.second_part;

        let mut rhs_is_y = rhs.pauli_letters.first_part.clone();
        *rhs_is_y &= &rhs.pauli_letters.second_part;
        let mut rhs_is_x = rhs.pauli_letters.second_part.negate();
        *rhs_is_x &= &rhs.pauli_letters.first_part;
        let mut rhs_is_z = rhs.pauli_letters.first_part.negate();
        *rhs_is_z &= &rhs.pauli_letters.second_part;

        let mut extra_factors_of_i: u8 = 0;

        {
            let mut self_is_x_rhs_is_y = self_is_x.clone();
            *self_is_x_rhs_is_y &= &rhs_is_y;
            let self_is_x_rhs_is_y = self_is_x_rhs_is_y.count_ones();
            extra_factors_of_i += u8::try_from(self_is_x_rhs_is_y % 4).unwrap();
        }
        {
            let mut self_is_x_rhs_is_z = self_is_x.clone();
            *self_is_x_rhs_is_z &= &rhs_is_z;
            let self_is_x_rhs_is_z = self_is_x_rhs_is_z.count_ones();
            extra_factors_of_i += u8::try_from((self_is_x_rhs_is_z % 4) * 3).unwrap();
        }

        {
            let mut self_is_y_rhs_is_z = self_is_y.clone();
            *self_is_y_rhs_is_z &= &rhs_is_z;
            let self_is_y_rhs_is_z = self_is_y_rhs_is_z.count_ones();
            extra_factors_of_i += u8::try_from(self_is_y_rhs_is_z % 4).unwrap();
        }
        {
            let mut self_is_y_rhs_is_x = self_is_y.clone();
            *self_is_y_rhs_is_x &= &rhs_is_x;
            let self_is_y_rhs_is_x = self_is_y_rhs_is_x.count_ones();
            extra_factors_of_i += u8::try_from((self_is_y_rhs_is_x % 4) * 3).unwrap();
        }

        {
            let mut self_is_z_rhs_is_x = self_is_z.clone();
            *self_is_z_rhs_is_x &= &rhs_is_x;
            let self_is_z_rhs_is_x = self_is_z_rhs_is_x.count_ones();
            extra_factors_of_i += u8::try_from(self_is_z_rhs_is_x % 4).unwrap();
        }
        {
            let mut self_is_z_rhs_is_y = self_is_z.clone();
            *self_is_z_rhs_is_y &= &rhs_is_y;
            let self_is_z_rhs_is_y = self_is_z_rhs_is_y.count_ones();
            extra_factors_of_i += u8::try_from((self_is_z_rhs_is_y % 4) * 3).unwrap();
        }

        self.overall_factor += extra_factors_of_i % 4;
        self.pauli_letters += &rhs.pauli_letters;
    }
}

pub struct SGenerators {
    generating_ms: Vec<PauliString>,
}

impl SGenerators {
    /// Create a new set of generators for an abelian subgroup of `P_n`
    /// Check that the generators actually all commute only if the `verify` flag says so.
    ///
    /// # Panics
    ///
    /// If we are doing verification and either there is not a consistent `n` for the number of qubits
    /// or they don't actually commute.
    ///
    pub fn new(generating_ms: Vec<PauliString>, verify: bool) -> Self {
        if verify {
            let mut num_qubits_each = generating_ms.iter().map(|mi| mi.pauli_letters.my_n);
            let all_should_be_this = num_qubits_each.next().unwrap();
            assert!(num_qubits_each.all(|num_qubits_now| num_qubits_now == all_should_be_this));
            for idx in 0..generating_ms.len() {
                for jdx in 0..idx {
                    assert!(
                        generating_ms[idx].commutes_with(&generating_ms[jdx]),
                        "{} does not commute with {}",
                        generating_ms[idx],
                        generating_ms[jdx]
                    );
                }
            }
        }
        Self { generating_ms }
    }

    pub fn num_qubits(&self) -> usize {
        self.generating_ms[0].pauli_letters.my_n
    }

    /// Assuming these were `n-k` independent generators
    ///
    /// # Panics
    ///
    /// If instead, we are verifying independence of the generators,
    /// then panic if they are not actually independent.
    pub fn my_k(&self, verify: bool) -> usize {
        if verify {
            let as_matrix = self.as_g_matrix();
            assert_eq!(as_matrix.rank(), self.generating_ms.len());
        }
        self.num_qubits() - self.generating_ms.len()
    }

    /// Is it giving a maximal isotropic subspace
    pub fn is_lagrangian(&self) -> bool {
        self.num_qubits() == self.generating_ms.len()
    }

    pub fn as_g_matrix(&self) -> BitMatrix {
        let num_qubits = self.num_qubits();
        BitMatrix::vstack_from_owned_iter(self.generating_ms.iter().map(|mi| {
            let mut left_part = BitMatrix::zeros(1, num_qubits);
            left_part.add_bits_to_row(&mi.pauli_letters.first_part, 0);
            let mut right_part = BitMatrix::zeros(1, num_qubits);
            right_part.add_bits_to_row(&mi.pauli_letters.second_part, 0);
            left_part.hstack(&right_part)
        }))
    }

    /// Changes the presentation of the group
    /// according to Gaussian elimination on the G-matrix.
    pub fn change_presentation(&mut self, blocksize: usize) {
        let mut as_matrix = self.as_g_matrix();
        as_matrix.gauss_with_proxy(true, blocksize, self);
        debug_assert_eq!(self.as_g_matrix(), as_matrix);
    }

    pub fn syndrome(&self, error: &PauliString) -> BitData {
        BitData::from(
            self.generating_ms
                .iter()
                .map(|mi| !mi.commutes_with(error))
                .collect::<Vec<_>>(),
        )
    }

    /// The weight 1 Pauli strings that centralize this subgroup
    /// They might be within this subgroup already or not
    #[allow(clippy::missing_panics_doc)]
    pub fn weight_one_centralizers(&self) -> impl Iterator<Item = PauliString> {
        let num_qubits = self.num_qubits();
        let mut to_return = Vec::with_capacity(num_qubits);
        let as_matrix = self.as_g_matrix();
        let num_rows = as_matrix.rows();
        let zero_col = BitMatrix::zeros(num_rows, 1);
        for idx in 0..num_qubits {
            let idxth_col = as_matrix.select_cols(&[idx]).expect("This column exists");
            let idx_plus_nth_col = as_matrix
                .select_cols(&[idx + num_qubits])
                .expect("This column exists");
            let idxth_is_zero = idxth_col == zero_col;
            let idx_plus_nth_is_zero = idx_plus_nth_col == zero_col;
            let both_cols_same = idxth_col == idx_plus_nth_col;
            if idxth_is_zero && idx_plus_nth_is_zero {
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::X,
                    IPower::One,
                ));
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::Y,
                    IPower::One,
                ));
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::Z,
                    IPower::One,
                ));
            } else if idxth_is_zero {
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::Z,
                    IPower::One,
                ));
            } else if idx_plus_nth_is_zero {
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::X,
                    IPower::One,
                ));
            } else if both_cols_same {
                to_return.push(PauliString::weight_one(
                    idx,
                    num_qubits,
                    PauliLetter::Y,
                    IPower::One,
                ));
            }
        }
        to_return.into_iter().flat_map(|with_no_prefactor| {
            let mut to_return = [
                with_no_prefactor.clone(),
                with_no_prefactor.clone(),
                with_no_prefactor.clone(),
                with_no_prefactor,
            ];
            to_return[1] *= IPower::PlusI;
            to_return[2] *= IPower::NegOne;
            to_return[3] *= IPower::NegI;
            to_return
        })
    }
}

impl RowOps for SGenerators {
    fn add_row(&mut self, from: usize, to: usize) {
        let from_piece = self.generating_ms[from].clone();
        self.generating_ms[to] *= &from_piece;
    }

    fn swap_rows(&mut self, from: usize, to: usize) {
        self.generating_ms.swap(from, to);
    }
}

mod test {

    #[test]
    fn pauli_1() {
        use super::{IPower, PauliString};
        let x_pauli = PauliString::parse("X").expect("Manifestly a Pauli String");
        let y_pauli = PauliString::parse("Y").expect("Manifestly a Pauli String");
        let z_pauli = PauliString::parse("Z").expect("Manifestly a Pauli String");
        let xyz_pauli = [x_pauli, y_pauli, z_pauli];
        for cur_squaring in &xyz_pauli {
            let mut cur_testing = cur_squaring.clone();
            cur_testing *= cur_squaring;
            assert_eq!(
                cur_testing,
                PauliString::parse("I").expect("Manifestly a Pauli String")
            );
        }
        for (which_left, cur_left) in xyz_pauli.iter().enumerate() {
            for (which_right, cur_right) in xyz_pauli.iter().enumerate() {
                if which_right == which_left {
                    continue;
                }
                let mut cur_testing = cur_left.clone();
                cur_testing *= cur_right;
                let which_expected = [0, 1, 2]
                    .into_iter()
                    .find(|idx| ![which_left, which_right].contains(idx))
                    .unwrap();
                if (which_right + 3 - which_left) % 3 == 1 {
                    cur_testing *= IPower::NegI;
                } else {
                    cur_testing *= IPower::PlusI;
                }
                assert_eq!(
                    cur_testing, xyz_pauli[which_expected],
                    "{}*{} vs {}",
                    cur_left, cur_right, xyz_pauli[which_expected]
                );
            }
        }
    }

    #[test]
    fn pauli_string_mul() {
        use super::PauliString;
        let mut m1 = PauliString::parse("XZZIX").expect("Manifestly a Pauli String");
        let m0 = PauliString::parse("IIIII").expect("Manifestly a Pauli String");
        m1 *= &m0;
        assert_eq!(
            m1,
            PauliString::parse("XZZIX").expect("Manifestly a Pauli String")
        );

        let m2 = PauliString::parse("ZXIZX").expect("Manifestly a Pauli String");
        m1 *= &m2;
        assert_eq!(
            m1,
            PauliString::parse("YYZZI").expect("Manifestly a Pauli String")
        );
    }

    #[test]
    fn stabilizer_5_1_3() {
        use super::PauliString;
        let m1 = PauliString::parse("XZZIX").expect("Manifestly a Pauli String");
        let m2 = PauliString::parse("ZXIZX").expect("Manifestly a Pauli String");
        let m3 = PauliString::parse("IZXZY").expect("Manifestly a Pauli String");
        let m4 = PauliString::parse("ZIZXY").expect("Manifestly a Pauli String");
        let ms = [m1, m2, m3, m4];
        let expected_weights = [4, 4, 4, 4];
        for (cur_m, expected_weight) in ms.iter().zip(expected_weights) {
            let as_symp = &cur_m.pauli_letters;
            assert_eq!(cur_m.weight(), expected_weight);
            assert_eq!(as_symp.weight(), expected_weight);
        }
        for idx in 0..4 {
            let ms_idx = &ms[idx];
            for jdx in 0..=idx {
                assert!(
                    ms_idx.commutes_with(&ms[jdx]),
                    "{} and {}",
                    ms_idx.pauli_letters,
                    ms[jdx].pauli_letters,
                );
            }
        }
    }

    #[test]
    fn shor_code() {
        use crate::{PauliString, SGenerators};
        let shor_shorthand = [
            "ZZI III III",
            "IZZ III III",
            "III ZZI III",
            "III IZZ III",
            "III III ZZI",
            "III III IZZ",
            "XXX XXX III",
            "III XXX XXX",
        ];
        let mut shor_code = SGenerators::new(
            shor_shorthand
                .into_iter()
                .map(|s| PauliString::parse(s).unwrap())
                .collect(),
            true,
        );

        let old_string = (0..8)
            .map(|idx| format!("{}\n", shor_code.generating_ms[idx]))
            .collect::<String>();
        shor_code.change_presentation(1);
        // can see that it has change the generators, but not the group
        let new_string = (0..8)
            .map(|idx| format!("{}\n", shor_code.generating_ms[idx]))
            .collect::<String>();
        assert_ne!(old_string, new_string);
    }
}
