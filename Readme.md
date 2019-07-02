[![pipeline status](https://vcs.eva.mpg.de/christian_heide/thrust/badges/master/pipeline.svg)](https://vcs.eva.mpg.de/christian_heide/thrust/commits/master) 
[![coverage report](https://vcs.eva.mpg.de/christian_heide/thrust/badges/master/coverage.svg)](https://vcs.eva.mpg.de/christian_heide/thrust/commits/master)

Apparently, GitLab-CI pipelines currently fail due to storage size reasons - not because of failing unit tests.

# thrust

This is another attempt to write a fast experimental ancient DNA damage aware short read mapper. The name "thrust" is 
an hommage to Udo Stenzel's "R-Candy", in which tradition it stands, and Rust, the programming language used for this 
project. This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

Thrust is based on the bidirectional FMD-index ([Li, 2012](https://academic.oup.com/bioinformatics/article/28/14/1838/218887)) 
with backtracking and lower-bound pruning of the search space. 
Improved algorithms and error models will be incorporated step by step as needed. 

Ancient DNA damage models can be included via the `SequenceDifferenceModel` trait. 

## Build and Install

Besides Rust, no additional dependencies are needed to compile. 

0. Install/update Rust (locally, to the current user's home directory):

    `curl https://sh.rustup.rs -sSf | sh`

1. Git clone thrust or simply download a release zip file:

    `git clone https://vcs.eva.mpg.de/christian_heide/thrust.git && cd thrust`

3. Build:

    `cargo build --release`
    
    The resulting binary file `thrust` is now in the subfolder `target/release/`.

4. Run!

    `cd target/release/`

    `./thrust index --reference /path/to/reference/hg19.fasta`
    
    `./thrust map --reads /path/to/reads/reads.fastq --reference /path/to/reference/hg19.fasta --output out.bam`

###### Optional
The replacement of step 3 with one of the following commands leads to increased performance on supported CPUs.

To build explicitly with SIMD support (should be available on most CPUs) use:

`cargo build --release --features simd-accel`

For AVX support (on recent CPUs like Intel Core i3/i5/i7 or recent AMD ones) use:

`RUSTFLAGS="-C target-cpu=native" cargo build --release --features "simd-accel avx-accel"`

To increase its verbosity, invoke the program like this:

`thrust -vvv index ...` or `thrust -vvv map ...`

## Performance/ Hardware Requirements

In #5 (tracking issue) benchmark results for mapping are shown.

Index generation for the human reference genome (hg19) unfortunately eats about 160GB of RAM (can certainly be improved 
by sampling the suffix array to k=32). Overall, the performance is comparable to `bwa aln` using ancient parameters. 
However, as soon as sampled suffix arrays are implemented, the performance will likely decrease a bit in favor of memory 
consumption.  

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
- [x] Incorporation of one or more of the ancient DNA damage models
  - [x] Add framework to specify PSSMs in the code
  - [x] simply penalize C-T changes according to Vindija pattern
  - [ ] Briggs et al. (2007)
  - [ ] ~~Skoglund et al. (2014)~~
  - [ ] Peyrégne (unpublished)
- [x] Recursive -> iterative k-mismatch-search
- [x] Extend match starting from the presumably lowest deaminated region (center) of a read
- [x] Calculate alignment score (with respect to damage/difference pattern, see above)
- [x] Move away from a stack-like data structure to a priority-queue (ordered by alignment scores) for partial matches
- [ ] Revisit mapping quality estimation (base-quality aware)
- [ ] BAM-IO
  - [x] BAM output
  - [ ] BAM input
- [x] Generate CIGAR string
- [ ] Evaluate precision, FPR, FNR
- [x] Only push to heap if necessary (perhaps check conditions in (inlined) function)
- [x] Problem: BFS does not work well with negative scores which we get with log-probs (-> slow). Perhaps add +1?
- [ ] Implement handy BWT construction algorithm on top of SA-IS (Okanohara & Sadakane, 2009)
- [ ] Derive sampled Suffix Arrays from BWT
- [x] Multithreading
- [ ] Cluster-enabled version
- [ ] _Paired-end sequencing_
