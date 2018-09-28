# Thrust
## We're in pre-alpha state. Basically nothing works at the moment, but there's hope...

This is an attempt to write an experimental ancient-DNA-damage-aware read-mapper. The name "Thrust" is an hommage to Udo Stenzel's "R-Candy", in which tradition it stands, and Rust, the programming language used for this project. This work depends heavily on the excellent [rust-bio](https://rust-bio.github.io/) crate. 

### TODO
- Use the FM*D*-index
- Search both strands
- Inexact matching
- Incorporation of an ancient DNA damage model
- Mapping quality estimation
- Multithreading
- BAM-IO
- Paired-end sequencing?
