use crate::{bitmatrix::*, bitvec::*};

use rayon::prelude::*;

pub trait ParallelMatrixOps {
    fn par_gauss(&mut self, full: bool);
}

impl ParallelMatrixOps for BitMatrix {
    fn par_gauss(&mut self, full: bool) {
        let chunk_size = 64 * self.col_blocks;
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
                }

                let row_vec = self.row(row).to_vec();

                let sz = self.data.len();
                self.data[((row1 + 1) * self.col_blocks)..sz]
                    .par_chunks_mut(chunk_size)
                    .for_each(|target_chunk| {
                        for row_start in (0..target_chunk.len()).step_by(self.col_blocks) {
                            if target_chunk.bit(BLOCKSIZE * row_start + pcol) {
                                target_chunk.xor_in(&row_vec, row_start);
                            }
                        }
                    });

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
                self.data[0..(row * self.col_blocks)]
                    .par_chunks_mut(chunk_size)
                    .for_each(|target_chunk| {
                        for row_start in (0..target_chunk.len()).step_by(self.col_blocks) {
                            if target_chunk.bit(BLOCKSIZE * row_start + pcol) {
                                target_chunk.xor_in(&row_vec, row_start);
                            }
                        }
                    });
            }
        }
    }
}
