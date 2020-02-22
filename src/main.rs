use clap::{crate_description, crate_version, value_t, App, AppSettings, Arg, SubCommand};

use mapad::{
    distributed::{dispatcher, worker},
    index, map,
    mismatch_bound::{Continuous, Discrete},
    sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
    utils,
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

    let matches = App::new(map::CRATE_NAME)
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
                        .takes_value(true)
                        .value_name("FASTA FILE"),
                )
                .arg(
                    Arg::with_name("seed")
                        .required(true)
                        .long("seed")
                        .help("Ambiguous bases are substituted with random ones. Two independent runs will lead to the same index when the same seed is used.")
                        .takes_value(true)
                        .default_value("1234")
                        .value_name("INT"),
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
                        .takes_value(true)
                        .value_name("BAM FILE"),
                )
                .arg(
                    Arg::with_name("reference")
                        .required(true)
                        .short("g")
                        .long("reference")
                        .help("Prefix of the file names of the index files. The reference FASTA file itself does not need to be present.")
                        .takes_value(true)
                        .value_name("FASTA FILE"),
                )
                .arg(
                    Arg::with_name("output")
                        .required(true)
                        .short("o")
                        .long("output")
                        .help("Path to output BAM file")
                        .takes_value(true)
                        .value_name("BAM FILE"),
                )
                .arg(
                    Arg::with_name("poisson_prob")
                        .short("p")
                        .group("allowed_mm")
                        .help("Minimum probability of the number of mismatches under 0.02 base error rate")
                        .takes_value(true)
                        // .default_value("0.04")
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                Arg::with_name("as_cutoff")
                    .short("c")
                    .group("allowed_mm")
                    .help("Per-base average log-likelihood cutoff")
                    .takes_value(true)
                    .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("as_cutoff_exponent")
                        .short("e")
                        .help("Exponent for read length dependency")
                        .takes_value(true)
                        .default_value("1.0")
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("library")
                        .required(true)
                        .short("l")
                        .long("library")
                        .help("Library preparation method")
                        .takes_value(true)
                        // .default_value("single_stranded")
                        .possible_values(&["single_stranded", "double_stranded"])
                        .value_name("METHOD")
                )
                .arg(
                    Arg::with_name("five_prime_overhang")
                        .required(true)
                        .short("f")
                        .help("5' overhang length parameter")
                        .takes_value(true)
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("three_prime_overhang")
                        .required_if("library", "single_stranded")
                        .short("t")
                        .help("3' overhang length parameter")
                        .takes_value(true)
                        .value_name("PROBABILITY")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ds_deamination_rate")
                        .required(true)
                        .short("d")
                        .help("Deamination rate in double-stranded stem of a read")
                        .takes_value(true)
                        // .default_value("0.02")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ss_deamination_rate")
                        .required(true)
                        .short("s")
                        .help("Deamination rate in single-stranded ends of a read")
                        .takes_value(true)
                        // .default_value("0.45")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("divergence")
                        .required(true)
                        .short("D")
                        .help("Divergence rate of the reference and target organisms")
                        .takes_value(true)
                        // .default_value("0.005")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("indel_rate")
                        .required(true)
                        .short("i")
                        .help("Expected rate of indels between reads and reference")
                        .takes_value(true)
                        .default_value("0.00001")
                        .value_name("RATE")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("chunk_size")
                        .required(true)
                        .long("batch_size")
                        .help("The number of reads that are processed in parallel")
                        .takes_value(true)
                        .default_value("250000")
                        .value_name("INT")
                )
                .arg(
                    Arg::with_name("dispatcher")
                        .long("dispatcher")
                        .help("Run in dispatcher mode for distributed computing in a network. Needs workers to be spawned externally to distribute work among them.")
                        .takes_value(false)
                )
                .arg(
                    Arg::with_name("port")
                        .required(true)
                        .long("port")
                        .help("Port to listen on in dispatcher mode")
                        .takes_value(true)
                        .default_value("3130"),
                ),

        )
        .subcommand(
            SubCommand::with_name("worker")
                .about("Spawns worker")
                .version(crate_version!())
                .arg(
                    Arg::with_name("host")
                        .required(true)
                        .long("host")
                        .help("Hostname or IP address of the running dispatcher node")
                        .takes_value(true)
                )
                .arg(
                    Arg::with_name("port")
                        .required(true)
                        .long("port")
                        .help("Port number of running dispatcher")
                        .takes_value(true)
                        .default_value("3130")
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
            if let Err(e) = index::run(
                index_matches.value_of("reference").unwrap(),
                value_t!(index_matches.value_of("seed"), u64).unwrap(),
            ) {
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

            let reads_path = map_matches.value_of("reads").unwrap();
            let reference_path = map_matches.value_of("reference").unwrap();
            let out_file_path = map_matches.value_of("output").unwrap();

            if map_matches.is_present("poisson_prob") {
                let mismatch_bound = Discrete::new(
                    value_t!(map_matches.value_of("poisson_prob"), f32)
                        .unwrap_or_else(|e| e.exit()),
                    0.02,
                    difference_model.get_representative_mismatch_penalty(),
                );

                let alignment_parameters = utils::AlignmentParameters {
                    penalty_gap_open: value_t!(map_matches.value_of("indel_rate"), f32)
                        .unwrap_or_else(|e| e.exit())
                        .log2(),
                    penalty_gap_extend: difference_model.get_representative_mismatch_penalty(), // FIXME
                    difference_model,
                    chunk_size: value_t!(map_matches.value_of("chunk_size"), usize)
                        .unwrap_or_else(|e| e.exit()),
                    mismatch_bound,
                };

                if map_matches.is_present("dispatcher") {
                    println!("Dispatcher mode");
                    let mut dispatcher = dispatcher::Dispatcher::new(
                        reads_path,
                        reference_path,
                        out_file_path,
                        &alignment_parameters,
                    )
                    .expect("Application error");

                    let port =
                        value_t!(map_matches.value_of("port"), u16).unwrap_or_else(|e| e.exit());
                    if let Err(e) = dispatcher.run(port) {
                        println!("Application error: {}", e);
                    }
                } else if let Err(e) = map::run(
                    reads_path,
                    reference_path,
                    out_file_path,
                    &alignment_parameters,
                ) {
                    println!("Application error: {}", e);
                }
            } else {
                let mismatch_bound = Continuous {
                    cutoff: value_t!(map_matches.value_of("as_cutoff"), f32)
                        .unwrap_or_else(|e| e.exit())
                        * -1.0,
                    exponent: value_t!(map_matches.value_of("as_cutoff_exponent"), f32).unwrap(),
                    representative_mismatch_penalty: difference_model
                        .get_representative_mismatch_penalty(),
                };

                let alignment_parameters = utils::AlignmentParameters {
                    penalty_gap_open: value_t!(map_matches.value_of("indel_rate"), f32)
                        .unwrap_or_else(|e| e.exit())
                        .log2(),
                    penalty_gap_extend: difference_model.get_representative_mismatch_penalty(), // FIXME
                    difference_model,
                    chunk_size: value_t!(map_matches.value_of("chunk_size"), usize)
                        .unwrap_or_else(|e| e.exit()),
                    mismatch_bound,
                };

                if map_matches.is_present("dispatcher") {
                    println!("Dispatcher mode");
                    let mut dispatcher = dispatcher::Dispatcher::new(
                        reads_path,
                        reference_path,
                        out_file_path,
                        &alignment_parameters,
                    )
                    .expect("Application error");

                    let port =
                        value_t!(map_matches.value_of("port"), u16).unwrap_or_else(|e| e.exit());
                    if let Err(e) = dispatcher.run(port) {
                        println!("Application error: {}", e);
                    }
                } else if let Err(e) = map::run(
                    reads_path,
                    reference_path,
                    out_file_path,
                    &alignment_parameters,
                ) {
                    println!("Application error: {}", e);
                }
            };
        }

        ("worker", Some(worker_matches)) => {
            let host = worker_matches.value_of("host").unwrap();
            let port = worker_matches.value_of("port").unwrap();
            let mut worker = worker::Worker::<SimpleAncientDnaModel, Discrete>::new(
                host,
                port,
            ).expect("Could not connect to dispatcher. Please check that it is running at the specified address.");

            worker.run().expect(
                "Could not successfully complete the tasks given. Please double-check the results.",
            );
        }
        _ => unreachable!(),
    }
}
