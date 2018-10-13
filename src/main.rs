extern crate clap;
#[macro_use]
extern crate log;
extern crate bio;
extern crate simple_logger;
extern crate thrust;

use bio::alphabets::Alphabet;
use clap::{App, AppSettings, Arg, SubCommand};

fn main() {
    let matches = App::new("Thrust")
        .about("An aDNA aware short-read mapper")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .arg(
            Arg::with_name("v")
                .short("v")
                .global(true)
                .multiple(true)
                .help("Sets the level of verbosity"),
        ).subcommand(
            SubCommand::with_name("index").arg(
                Arg::with_name("reference")
                    .required(true)
                    .long("reference")
                    .help("FASTA file containing the genome to be indexed")
                    .value_name("FASTA file"),
            ),
        ).subcommand(
            SubCommand::with_name("map")
//                .arg(
//                    Arg::with_name("index")
//                        .required(true)
//                        .long("index")
//                        .help("idx file of the genome we are about to map against")
//                        .value_name("idx file"),
//                ).arg(
//                    Arg::with_name("outfile")
//                        .required(true)
//                        .long("out")
//                        .help("output BAM file for the aligned reads")
//                        .value_name("BAM file"),
//                )
                .arg(
                    Arg::with_name("reads")
                        .required(true)
                        .long("reads")
                        .help("FASTQ file containing adapter-trimmed and quality-controlled reads")
                        .value_name("FASTQ file"),
                ),
        ).get_matches();

    simple_logger::init_with_level(match matches.occurrences_of("v") {
        0 => log::Level::Warn,
        1 => log::Level::Info,
        2 => log::Level::Debug,
        3 | _ => log::Level::Trace,
    }).unwrap();

    debug!("Rank-transform alphabet");
    let symbols = b"$ACGTN";
    let rank_symbols = symbols
        .iter()
        .enumerate()
        .map(|(i, _v)| i as u8)
        .collect::<Vec<_>>();

    // Define two alphabets, one with actual chars and
    // another one containing its corresponding ranks
    let bwt_alphabet = Alphabet::new(symbols.iter());
    let rank_alphabet = Alphabet::new(&rank_symbols);

    match matches.subcommand() {
        ("index", Some(index_matches)) => {
            if let Err(e) = thrust::index::run(
                index_matches.value_of("reference").unwrap(),
                &bwt_alphabet,
                &rank_alphabet,
            ) {
                println!("Application error: {}", e);
            }
        }
        ("map", Some(map_matches)) => {
            if let Err(e) = thrust::map::run(map_matches.value_of("reads").unwrap(), &bwt_alphabet)
            {
                println!("Application error: {}", e);
            }
        }
        _ => unreachable!(),
    }
}
