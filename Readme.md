[![pipeline status](https://vcs.eva.mpg.de/christian_heide/mapAD/badges/master/pipeline.svg)](https://vcs.eva.mpg.de/christian_heide/mapAD/commits/master) 
[![coverage report](https://vcs.eva.mpg.de/christian_heide/mapAD/badges/master/coverage.svg)](https://vcs.eva.mpg.de/christian_heide/mapAD/commits/master)

Apparently, GitLab-CI pipelines currently fail due to storage size reasons - not because of failing unit tests.

# mapAD <img src="https://vcs.eva.mpg.de/uploads/-/system/project/avatar/1089/480px-Uracil-3D-balls.png" alt="drawing" width="25"/>

This is another attempt to write a fast experimental ancient DNA damage aware short read mapper. 
This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

mapAD is based on the bidirectional FMD-index ([Li, 2012](https://academic.oup.com/bioinformatics/article/28/14/1838/218887)) 
with backtracking and lower-bound pruning of the search space. 
Improved algorithms and error models will be incorporated step by step as needed. 

Ancient DNA damage models can be included via the `SequenceDifferenceModel` trait. 
The default (and only) impl is based on Udo Stenzel's ANFO/r-candy. 

mapAD should be ready to test. Please note that it only reports one (the best) alignment per read so far.  
The program needs ~ 160GB RAM to align against the human reference genome hg19. The mapping speed is 
roughly comparable to BWA aln (ancient parameters).  
To parallelize alignments, it will use all idle CPU cores on the machine it is run on. 

## Build and Install

Besides Rust, no additional dependencies are needed to compile. 

0. Install/update Rust (locally, to the current user's home directory):

    `curl https://sh.rustup.rs -sSf | sh`

1. Git clone mapAD or simply download a release zip file:

    `git clone https://vcs.eva.mpg.de/christian_heide/mapAD.git && cd mapAD`

3. Build:

    `cargo build --release`
    
    The resulting binary file `mapad` is now in the subfolder `target/release/`.

4. Run!

    `cd target/release/`

    `./mapad index --reference /path/to/reference/hg19.fasta`
    
###### Optional
- The binary can be tailored to modern CPU architectures to increase its performance:

`RUSTFLAGS="-C target-cpu=haswell" cargo build --release"`

or

`RUSTFLAGS="-C target-cpu=native" cargo build --release"`

Please note that the resulting binary will probably only run on systems with the same or a newer CPU. 

- The number of `v`s passed to the program determines the level of verbosity:

`mapad -vvv index ...` or `mapad -vvv map ...`

## Usage
The subprograms `mapad index` and `mapad map` will index the reference and map reads to it, respectively. 
Adding ` --help` will print a list of available and required command line options. 

### Example
The following example aligns reads that are expected to have a Vindija-like deamination pattern to an existing index of the hg19 reference.
#### Local 
```bash
mapad -vvv map \
--library single_stranded \
-p 0.02 \
-f 0.475 \
-t 0.475 \
-d 0.001 \
-s 0.9 \
-D 0.02 \
-i 0.00001 \
--reads "${input_bam}" \
--reference "/mnt/scratch/chris/hg19_evan/whole_genome.fa" \
--output "${output_bam}"
```

#### Distributed
The following example starts a dispatcher node and then spawns multi-threaded workers on cluster nodes that have more than 64GB of free RAM. 
Start the dispatcher:
```bash
mapad -vv map \
--dispatcher \
# ... (see local example)
```
Spawn workers:
```bash
qsub -N "mapAD_worker" -pe "smp" 1-128 -t 1-65535 -l "h_vmem=64G,s_vmem=64G,virtual_free=64G,mem_free=64G,class=*" -j "y" -R "y" "mapad_worker.sh"
```

`mapad_worker.sh`:
```bash
#!/bin/bash
mapad -vvv worker --host <HOST>
```

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
  - [x] Model inspired by Udo Stenzel's ANFO/r-candy
  - [ ] Peyrégne (unpublished)
- [x] Recursive -> iterative k-mismatch-search
- [x] Extend match starting from the presumably lowest deaminated region (center) of a read
- [x] Calculate alignment score (with respect to damage/difference pattern, see above)
- [x] Move away from a stack-like data structure to a priority-queue (ordered by alignment scores) for partial matches
- [x] Revisit mapping quality estimation
- [x] BAM-IO
  - [x] BAM output
  - [x] BAM input
- [x] Generate CIGAR string
- [ ] Evaluate precision, FPR, FNR
- [x] Only push to heap if necessary (perhaps check conditions in (inlined) function)
- [x] Problem: BFS does not work well with negative scores which we get with log-probs (-> slow). Perhaps add +1?
- [ ] Implement handy BWT construction algorithm on top of SA-IS (Okanohara & Sadakane, 2009)
- [ ] Derive sampled Suffix Arrays from BWT
- [x] Multithreading
- [ ] Cluster-enabled version
- [ ] _Paired-end sequencing_
