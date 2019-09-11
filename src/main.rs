use clap::{crate_description, crate_version, value_t, App, AppSettings, Arg, SubCommand};

use mapad::{
    sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
    {index, map, utils},
};

fn main() {
    let probability_validator = |v: String| {
        let error_message = String::from("Please specify a value between 0 and 1");
        let v: f32 = match v.parse() {
            Ok(s) => s,
            Err(_) => return Err(error_message),
        };
        if (v >= 0.0) && (v <= 1.0) {
            return Ok(());
        }
        Err(error_message)
    };

    let matches = App::new("mapAD")
        .about(crate_description!())
        .version(crate_version!())
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .arg(
            Arg::with_name("v")
                .short("v")
                .global(true)
                .multiple(true)
                .help("Sets the level of verbosity"),
        )
        .subcommand(
            SubCommand::with_name("index")
                .about("Indexes a genome file")
                .version(crate_version!())
                .arg(
                    Arg::with_name("reference")
                        .required(true)
                        .short("g")
                        .long("reference")
                        .help("FASTA file containing the genome to be indexed")
                        .value_name("FASTA FILE"),
                ),
        )
        .subcommand(
            SubCommand::with_name("map")
                .about("Maps reads to an indexed genome")
                .version(crate_version!())
                .arg(
                    Arg::with_name("reads")
                        .required(true)
                        .short("r")
                        .long("reads")
                        .help("BAM file containing adapter-trimmed and quality-controlled reads")
                        .value_name("BAM FILE"),
                )
                .arg(
                    Arg::with_name("reference")
                        .required(true)
                        .short("g")
                        .long("reference")
                        .help("Prefix of the file names of the index files. The reference FASTA file itself does not need to be present.")
                        .value_name("FASTA FILE"),
                )
                .arg(
                    Arg::with_name("output")
                        .required(true)
                        .short("o")
                        .long("output")
                        .help("Path to output BAM file")
                        .value_name("BAM FILE"),
                )
                .arg(
                    Arg::with_name("poisson_prob")
                        .required(true)
                        .short("p")
                        .conflicts_with("max_diff")
                        // .default_value("0.04")
                        .help("Minimum probability of the number of mismatches under 0.02 base error rate")
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("library")
                        .required(true)
                        .short("l")
                        .long("library")
                        // .default_value("single_stranded")
                        .possible_values(&["single_stranded", "double_stranded"])
                        .help("Library preparation method")
                        .value_name("METHOD")
                )
                .arg(
                    Arg::with_name("five_prime_overhang")
                        .required(true)
                        .short("f")
                        .help("5' overhang length parameter")
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("three_prime_overhang")
                        .required_if("library", "single_stranded")
                        .short("t")
                        .help("3' overhang length parameter")
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ds_deamination_rate")
                        .required(true)
                        .short("d")
                        // .default_value("0.02")
                        .help("Deamination rate in double-stranded stem of a read")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ss_deamination_rate")
                        .required(true)
                        .short("s")
                        // .default_value("0.45")
                        .help("Deamination rate in single-stranded ends of a read")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("divergence")
                        .required(true)
                        .short("D")
                        // .default_value("0.005")
                        .help("Divergence rate of the reference and target organisms")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("indel_rate")
                        .required(true)
                        .short("i")
                        .default_value("0.00001")
                        .help("Expected rate of indels between reads and reference")
                        .value_name("RATE")
                        .validator(probability_validator),
                ),
        )
        .get_matches();

    simple_logger::init_with_level(match matches.occurrences_of("v") {
        0 => log::Level::Warn,
        1 => log::Level::Info,
        2 => log::Level::Debug,
        3 | _ => log::Level::Trace,
    })
    .unwrap();

    match matches.subcommand() {
        ("index", Some(index_matches)) => {
            if let Err(e) = index::run(index_matches.value_of("reference").unwrap()) {
                println!("Application error: {}", e);
            }
        }
        ("map", Some(map_matches)) => {
            let library_prep = match map_matches.value_of("library").unwrap() {
                "single_stranded" => LibraryPrep::SingleStranded {
                    five_prime_overhang: value_t!(map_matches.value_of("five_prime_overhang"), f32)
                        .unwrap_or_else(|e| e.exit()),
                    three_prime_overhang: value_t!(
                        map_matches.value_of("three_prime_overhang"),
                        f32
                    )
                    .unwrap_or_else(|e| e.exit()),
                },
                "double_stranded" => LibraryPrep::DoubleStranded(
                    value_t!(map_matches.value_of("five_prime_overhang"), f32)
                        .unwrap_or_else(|e| e.exit()),
                ),
                _ => unreachable!(),
            };
            let difference_model = SimpleAncientDnaModel {
                library_prep,
                ds_deamination_rate: value_t!(map_matches.value_of("ds_deamination_rate"), f32)
                    .unwrap_or_else(|e| e.exit()),
                ss_deamination_rate: value_t!(map_matches.value_of("ss_deamination_rate"), f32)
                    .unwrap_or_else(|e| e.exit()),
                // Divergence is divided by three because it is used for testing each of the three possible substitutions
                divergence: value_t!(map_matches.value_of("divergence"), f32)
                    .unwrap_or_else(|e| e.exit())
                    / 3.0,
            };
            let alignment_parameters = utils::AlignmentParameters {
                base_error_rate: 0.02,
                poisson_threshold: value_t!(map_matches.value_of("poisson_prob"), f64)
                    .unwrap_or_else(|e| e.exit()),
                penalty_gap_open: value_t!(map_matches.value_of("indel_rate"), f32)
                    .unwrap_or_else(|e| e.exit())
                    .log2(),
                penalty_gap_extend: difference_model.get_representative_mismatch_penalty(), // TODO
                difference_model,
            };
            if let Err(e) = map::run(
                map_matches.value_of("reads").unwrap(),
                map_matches.value_of("reference").unwrap(),
                map_matches.value_of("output").unwrap(),
                &alignment_parameters,
            ) {
                println!("Application error: {}", e);
            }
        }
        _ => unreachable!(),
    }
}
