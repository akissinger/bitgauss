[![pypi][]](https://pypi.org/project/bitgauss/)
[![py-version][]](https://pypi.org/project/bitgauss/)
[![crates][]](https://crates.io/crates/bitgauss)
[![docs][]](https://bitgauss.readthedocs.io/en/latest/)
[![rs-docs][]](https://docs.rs/bitgauss)

  [pypi]: https://img.shields.io/pypi/v/bitgauss
  [py-version]: https://img.shields.io/pypi/pyversions/bitgauss
  [crates]: https://img.shields.io/crates/v/bitgauss
  [docs]: https://app.readthedocs.org/projects/bitgauss/badge/?version=latest
  [rs-docs]: https://img.shields.io/docsrs/bitgauss?label=rust%20docs


bitgauss is a library for doing linear algebra over the 2-element finite field. It consists of a [Rust crate](https://docs.rs/bitgauss) with [Python bindings](https://bitgauss.readthedocs.io/en/latest/). To use the Rust crate, add `bitgauss` as a dependency to `Cargo.toml`, or to use the Python bindings, run:

```bash
pip install bitgauss
```

It provides a `BitMatrix` class, which is a bitpacked 2D matrix that implements fast linear algebraic operations using fast, vectorized bitwise operations. The main features are:
- getting and setting individual matrix elements (as `bool`s)
- fast row operations and dot product using bitwise operations
- fast in-place and out-of-place matrix transpose using a [recursive block method](https://github.com/dsnet/matrix-transpose)
- horizontal and vertical concatenation of matrices
- matrix multiplication
- Gaussian elimination and related methods (e.g. rank, inverse, and nullspace)

The goal of this library is to keep things small and simple, but if you need some core functionality that is missing, feel free to send a pull request!

Here is a simple example of using the `BitMatrix` class from Python:

```python
from bitgauss import BitMatrix

# Construct a 300x400 matrix whose entries are given by the bool-valued function
m1 = BitMatrix.build(300, 400, lambda i, j: (i + j) % 2 == 0)

# Construct a random 80x300 matrix with an optional random seed
m2 = BitMatrix.random(80, 300, seed=1)

# Construct a random invertible 300x300 matrix
m3 = BitMatrix.random_invertible(300, seed=1)

m4 = m2 * m3           # Matrix multiplication
m3_inv = m3.inverse()  # Returns the inverse
m1_t = m1.transposed() # Returns transpose
m1.transpose_inplace() # Transpose inplace (padding if necessary)
m1.gauss()             # Transform to row-echelon form
m1.gauss(full=True)    # Transform to reduced row-echelon form
ns = m1.nullspace()    # Returns a spanning set for the nullspace
```

...or from Rust:

```rust
use bitmatrix::BitMatrix;
use std::ops::Mul;
use rand::{Rng, SeedableRng, rngs::SmallRng};

let mut rng = SmallRng::seed_from_u64(1);

// Construct a 300x400 matrix whose entries are given by the bool-valued function
let mut m1 = BitMatrix::build(300, 400, |i, j| (i + j) % 2 == 0);

// Construct a random 80x300 matrix using the given random number generator
let m2 = BitMatrix::random(&mut rng, 80, 300);

// Construct a random invertible 300x300 matrix
let m3 = BitMatrix::random_invertible(&mut rng, 300);

let m4 = &m2 * &m3;         // Matrix multiplication
let m3_inv = m3.inverse();  // Returns the inverse
let m1_t = m1.transposed(); // Returns transpose
m1.transpose_inplace();     // Transpose inplace (padding if necessary)
m1.gauss(false);            // Transform to row-echelon form
m1.gauss(true);             // Transform to reduced row-echelon form
let ns = m1.nullspace();    // Returns a spanning set for the nullspace
```

For more info, have a look at the [Python documentation](https://bitgauss.readthedocs.io/en/latest/), the [Rust documentation](https://docs.rs/bitgauss), or the Jupyter notebook [getting_started.ipynb](https://github.com/akissinger/bitgauss/blob/main/demo/getting_started.ipynb).

# Performance

The `BitMatrix` class is substantially faster than a pure Python implementation of a binary matrix without bitpacking. For example, in [getting_started.ipynb](https://github.com/akissinger/bitgauss/blob/main/demo/getting_started.ipynb) there is a calculation that is more than 100X faster using `bitgauss`:

```python
m = BitMatrix.random(1000, 2000, seed=42) # bitgauss
m1 = Mat2(m.to_int_list()) # pure python

%timeit m.copy().gauss()
# 82.2 ms ± 1.7 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

%timeit m1.copy().gauss()
# 10.1 s ± 172 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

The difference is probably less pronounced with `numpy` or other fast compiled linear algebra backends. If you'd like to contribute some benchmarks vs. those, feel free.

## A note on SIMD

`bitgauss` is designed to take advantage of the fact that modern hardware and compilers are able to process many bit-entries at once using [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data). While this package doesn't explicitly use SIMD (via e.g. the experimental [std::simd](https://doc.rust-lang.org/std/simd/index.html) module in nightly Rust), it relies on the fact that [LLVM](https://llvm.org/), which the Rust compiler uses under the hood, is already very good at using SIMD automatically when needed and available. My experiments and benchmarking, this seems to almost always do better than introducing SIMD instructions by hand. Please keep this in mind before you spend a lot of time writing a PR that adds explicit SIMD instructions. Basically, don't do what I did. :)