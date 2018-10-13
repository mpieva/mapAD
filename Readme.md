# Thrust

_We're in pre-alpha state. Very few things work at the moment, but there's hope..._

This is another attempt to write a _fast_ experimental ancient DNA damage aware short read mapper. The name "Thrust" is 
an hommage to Udo Stenzel's "R-Candy", in which tradition it stands, and Rust, the programming language used for this 
project. This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

The first versions will basically be re-implementations of BWA (FM-index with backtracking and lower-bound pruning of 
the search space). Improved algorithms and error models will be incorporated step by step as needed. 

As of now, it's not quite clear whether or not to prefer a rather simple approach (as used in bisulfite sequencing) to 
the ancient DNA damage problem over a more complex error model.

### Performance/ Hardware Requirements
First tests suggest that index generation for the human reference genome (hg19) unfortunately eats about 12G of RAM 
(can certainly be improved,  but it's not an high-priority issue). 
Holding that index in memory during the mapping works fine on a 8G laptop. Indexing hg19 has a runtime of about 
20 minutes on fast machines. Most of the time is spent on writing/reading and (de-)serializing of index files.  

### To do
- [x] Use FM**D**-index
- [x] Search both strands
- [x] Introduce subcommands {index, map}
- [x] Save reference index to disk
- [x] Refactor crate/mod structure
- [ ] Make SIMD vectorization easily configurable at compile time (as "feature")
- [ ] Inexact matching
- [ ] Incorporation of one or more of the ancient DNA damage models 
  - [ ] simply ignore C-T deaminations
  - [ ] Briggs et al. (2007)
  - [ ] Skoglund et al. (2014)
  - [ ] Peyrégne (unpublished)
- [ ] Mapping quality estimation
- [ ] Multithreading
- [ ] BAM-IO
- [ ] _Cluster-enabled version_
- [ ] _Paired-end sequencing_
