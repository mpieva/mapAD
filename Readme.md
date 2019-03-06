[![pipeline status](https://vcs.eva.mpg.de/christian_heide/thrust/badges/master/pipeline.svg)](https://vcs.eva.mpg.de/christian_heide/thrust/commits/master) 
[![coverage report](https://vcs.eva.mpg.de/christian_heide/thrust/badges/master/coverage.svg)](https://vcs.eva.mpg.de/christian_heide/thrust/commits/master)

# thrust

This is another attempt to write a fast experimental ancient DNA damage aware short read mapper. The name "thrust" is 
an hommage to Udo Stenzel's "R-Candy", in which tradition it stands, and Rust, the programming language used for this 
project. This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

The first versions will basically be re-implementations of BWA (FM-index with backtracking and lower-bound pruning of 
the search space). Improved algorithms and error models will be incorporated step by step as needed. 

As of now, it's not quite clear whether or not to prefer a rather simple approach (as used in bisulfite sequencing) to 
the ancient DNA damage problem over a more complex error model.

## Build and Install

Besides Rust, no additional dependencies are needed to compile. 

0. Install/update Rust (locally, to the current user's home directory):

    `curl https://sh.rustup.rs -sSf | sh`

1. Git clone thrust or simply download a release zip file:

    `git clone https://vcs.eva.mpg.de/christian_heide/thrust.git && cd thrust`

3. Build:

    `cargo build --release`
    
    The resulting binary file `thrust` is now in the subfolder `./target/release/`.

4. Run!

    `./target/release/thrust index --reference ./example/chr22.fa`
    
    `./target/release/thrust map --reads ./example/simulated_reads/test.bwa.1.fa`

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
- [x] Don't rank-transform during early development. Go back to upstream rust-bio
- [x] Inexact matching
- [x] It's time to start testing!
- [ ] Incorporation of one or more of the ancient DNA damage models
  - [x] simply ignore C-T deaminations
  - [ ] simply penalize C-T changes according to Vindija pattern
  - [ ] Briggs et al. (2007)
  - [ ] Skoglund et al. (2014)
  - [ ] Peyrégne (unpublished)
- [x] Recursive -> iterative k-mismatch-search
- [x] Extend match starting from the presumably lowest deaminated region (center) of a read
- [ ] Mapping quality estimation
- [ ] BAM-IO
- [ ] Implement handy BWT construction algorithm on top of SA-IS (Okanohara & Sadakane, 2009)
- [ ] Derive sampled Suffix Arrays from BWT
- [ ] Multithreading
- [ ] _Cluster-enabled version_
- [ ] _Paired-end sequencing_

## Performance/ Hardware Requirements

First tests suggest that index generation for the human reference genome (hg19) unfortunately eats about 160GB of RAM 
(can certainly be improved by sampling the suffix array to k=32). As long as rust-bio does not have a working suffix 
array sampling implementation, whole-genome mapping is utopic. As of now, tests are performed on chr22. 
