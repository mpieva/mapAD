# Thrust

_We're in pre-alpha state. Very few things work at the moment, but there's hope..._

This is another attempt to write a _fast_ experimental ancient DNA damage aware short read mapper. The name "Thrust" is 
an hommage to Udo Stenzel's "R-Candy", in which tradition it stands, and Rust, the programming language used for this 
project. This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

The first versions will basically be (slightly improved) re-implementations of BWA (FM-index with backtracking and lower-bound pruning of 
the search space). Improved algorithms and error models will be incorporated step by step as needed. 

As of now, it's not quite clear whether or not to prefer a rather simple approach (as used in bisulfite sequencing) to 
the ancient DNA damage problem over a more complex error model.

## Build and Install

Besides Rust, no additional dependencies are needed to compile. 
Rust must be installed to build from source, see: https://www.rust-lang.org/.

1. Git clone or simply download a release zip file:

    `git clone repository/thrust.git`

2. Enter the downloaded directory:

    `cd thrust`

3. Build:

    `cargo build --release`

###### Optional
The replacement of step 3 with one of the following commands leads to increased performance on supported CPUs.

To build explicitly with SIMD support (should be available on most CPUs) use:

`cargo build --release --features simd-accel`

for AVX support (on recent CPUs like Intel Core i3/i5/i7 or recent AMD ones) use:

`RUSTFLAGS="-C target-cpu=native" cargo build --release --features "simd-accel avx-accel"`

## To do

- [x] Use FM**D**-index
- [x] Search both strands
- [x] Introduce subcommands {index, map}
- [x] Save reference index to disk
- [x] Refactor crate/mod structure
- [x] Make SIMD vectorization easily configurable at compile time (as "feature")
- [x] Compress index files for hopefully faster IO (libflate)
- [ ] Inexact matching
- [ ] Incorporation of one or more of the ancient DNA damage models 
  - [ ] simply ignore C-T deaminations
  - [ ] Briggs et al. (2007)
  - [ ] Skoglund et al. (2014)
  - [ ] Peyrégne (unpublished)
- [ ] Mapping quality estimation
- [ ] Multithreading
- [ ] BAM-IO
- [ ] Fix suffix array sampling in rust-bio
- [ ] _Cluster-enabled version_
- [ ] _Paired-end sequencing_

## Performance/ Hardware Requirements

First tests suggest that index generation for the human reference genome (hg19) unfortunately eats about 160GB of RAM 
(can certainly be improved by sampling the suffix array to k=32). As long as rust-bio does not have a working suffix 
array sampling implementation, whole-genome mapping is utopic. As of now, tests are performed on chr22. 
