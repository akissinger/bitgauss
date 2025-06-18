bitgauss is a Rust library for doing linear algebra over the 2-element finite field. It's main data structure is `BitMatrix`, which a bitpacked 2D matrix that implements fast linear algebraic operations using bitwise operations.

Currently, its features are very minimal, but I expect these will grow over time. The main ones are:
- getting and setting individual matrix elements (as `bool`s)
- fast row operations and dot product using bitwise operations
- fast in-place and out-of-place matrix transpose using a [recursive block method](https://github.com/dsnet/matrix-transpose)
- horizontal and vertical concatenation of matrices
- matrix multiplication
- Gaussian elimination and related methods (e.g. rank and inverse)

An additional planned feature is computing spanning sets for images and kernels of matrices. Other features may be added as and when they are needed by other projects. If you need something, feel free to send a pull request!

Much of the speed comes from the fact that modern hardware and compilers are able to process many bit-entries at once using [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data). While this package doesn't explicitly use SIMD (via e.g. the experimental [std::simd](https://doc.rust-lang.org/std/simd/index.html) module in nightly Rust), it relies on the fact that [LLVM](https://llvm.org/), which the Rust compiler uses under the hood, is already very good at using SIMD automatically when needed and available. From some basic benchmarking, this seems to always do a bit better than introducing SIMD instructions by hand.