bitgauss is a Rust library for doing linear algebra over the 2-element finite field. It's main data structure is `BitMatrix`, which a bitpacked 2D matrix that implements fast linear algebraic operations using bitwise operations.

Currently, its features are very minimal, but I expect these will grow over time. The main ones are:
- getting and setting individual matrix elements (as `bool`s)
- fast row operations and dot product using bitwise operations
- fast in-place and out-of-place matrix transpose using a [recursive block method](https://github.com/dsnet/matrix-transpose)
- matrix multiplication
- Gaussian elimination and related methods (e.g. rank and inverse)

Some planned features are horizontal/vertical stacking and splitting and computing spanning sets for images and kernels of matrices. Other features may be added as and when they are needed by other projects. If you need something, feel free to send a pull request!

Much of the speed comes from the fact that modern hardware and compilers are able to process many bit-entries at once using [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data). I have done some experiments with introducing SIMD operations manually into this code using the experimental [std::simd](https://doc.rust-lang.org/std/simd/index.html) module in nightly Rust. However, [LLVM](https://llvm.org/), which the Rust compiler uses under the hood, is already very good at using SIMD automatically when needed and available, and I found that it was nearly always fast just to let LLVM do its magic rather than invoking SIMD instructions by hand.