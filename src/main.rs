use clap::{
    crate_description, crate_version, value_t, App, AppSettings, Arg, ArgMatches, SubCommand,
};
use log::{error, info, warn};
use simple_logger::SimpleLogger;

use mapad::{
    distributed::{dispatcher, worker},
    index, map,
    mismatch_bounds::{Continuous, Discrete},
    sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
    utils::AlignmentParameters,
};

fn main() {
    handle_arguments(define_cli());
}

fn define_cli<'a>() -> ArgMatches<'a> {
    let probability_validator = |v: String| {
        let error_message = String::from("Please specify a value between 0 and 1");
        let v: f32 = match v.parse() {
            Ok(s) => s,
            Err(_) => return Err(error_message),
        };
        if (0.0..=1.0).contains(&v) {
            Ok(())
        } else {
            Err(error_message)
        }
    };

    App::new(map::CRATE_NAME)
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
                        .help("BAM or FASTQ file that contains the reads to be aligned")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("reference")
                        .required(true)
                        .short("g")
                        .long("reference")
                        .help("Prefix of the file names of the index files. The reference FASTA file itself does not need to be present.")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output")
                        .required(true)
                        .short("o")
                        .long("output")
                        .help("Path to output BAM file")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("poisson_prob")
                        .short("p")
                        .group("allowed_mm")
                        .help("Minimum probability of the number of mismatches under `-D` base error rate")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("as_cutoff")
                        .short("c")
                        .group("allowed_mm")
                        .help("Per-base average alignment score cutoff (-c > AS / read_len^e ?)")
                        .takes_value(true)
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("as_cutoff_exponent")
                        .short("e")
                        .help("Exponent to be applied to the read length (ignored if `-c` is not used)")
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
                        .possible_values(&["single_stranded", "double_stranded"])
                        .value_name("STRING")
                )
                .arg(
                    Arg::with_name("five_prime_overhang")
                        .required(true)
                        .short("f")
                        .help("5'-overhang length parameter")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("three_prime_overhang")
                        .required_if("library", "single_stranded")
                        .short("t")
                        .help("3'-overhang length parameter")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ds_deamination_rate")
                        .required(true)
                        .short("d")
                        .help("Deamination rate in double-stranded stem of a read")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("ss_deamination_rate")
                        .required(true)
                        .short("s")
                        .help("Deamination rate in single-stranded ends of a read")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("divergence")
                        .required(true)
                        .short("D")
                        .help("Divergence / base error rate")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .default_value("0.02")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("indel_rate")
                        .required(true)
                        .short("i")
                        .help("Expected rate of indels between reads and reference")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::with_name("gap_extension_penalty")
                        .required(true)
                        .short("x")
                        .help("Gap extension penalty as a fraction of the representative mismatch penalty")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .default_value("1.0")
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
                    Arg::with_name("ignore_base_quality")
                        .long("ignore_base_quality")
                        .help("Ignore base qualities in scoring models")
                        .takes_value(false)
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
                        .default_value("3130")
                        .value_name("INT")
                )
                .arg(
                    Arg::with_name("num_threads")
                        .required(true)
                        .long("threads")
                        .help(&format!("Maximum number of threads. If 0 or unspecified, {} will select the number of threads automatically.", map::CRATE_NAME))
                        .takes_value(true)
                        .default_value("0")
                        .value_name("INT")
                )
                .arg(
                    Arg::with_name("gap_dist_ends")
                        .required(true)
                        .long("gap_dist_ends")
                        .help("Disallow allow gaps at read ends (configurable range)")
                        .takes_value(true)
                        .default_value("5")
                        .value_name("INT")
                )
                .arg(
                    Arg::with_name("stack_limit_abort")
                        .long("stack_limit_abort")
                        .help("Abort alignment when stack size limit is reached instead of trying to recover.")
                        .takes_value(false)
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
                        .value_name("HOST")
                )
                .arg(
                    Arg::with_name("port")
                        .required(true)
                        .long("port")
                        .help("Port number of running dispatcher")
                        .takes_value(true)
                        .default_value("3130")
                        .value_name("INT")
                )
                .arg(
                    Arg::with_name("num_threads")
                        .required(true)
                        .long("num_threads")
                        .help(&format!("Maximum number of threads. If 0 or unspecified, {} will select the number of threads automatically.", map::CRATE_NAME))
                        .takes_value(true)
                        .default_value("0")
                        .value_name("INT")
                ),
        )
        .get_matches()
}

fn handle_arguments(matches: ArgMatches) {
    SimpleLogger::new()
        .with_level(match matches.occurrences_of("v") {
            0 => log::LevelFilter::Info,
            1 => log::LevelFilter::Debug,
            _ => log::LevelFilter::Trace,
        })
        .with_utc_timestamps()
        .init()
        .expect("This is not expected to fail");
    warn!("Use UTC timestamps for log entries");

    match matches.subcommand() {
        ("index", Some(arg_matches)) => {
            start_indexer(arg_matches);
        }
        ("map", Some(arg_matches)) => {
            start_mapper(arg_matches);
        }
        ("worker", Some(arg_matches)) => {
            start_worker(arg_matches);
        }
        _ => unreachable!(),
    }
}

fn start_indexer(arg_matches: &ArgMatches) {
    let seed = value_t!(arg_matches.value_of("seed"), u64).unwrap_or_else(|e| e.exit());
    let reference_path = arg_matches
        .value_of("reference")
        .expect("Presence is ensured by CLI definition");

    if let Err(e) = index::run(reference_path, seed) {
        error!("{}", e);
    }
}

fn start_mapper(map_matches: &ArgMatches) {
    let reads_path = map_matches
        .value_of("reads")
        .expect("Presence is ensured by CLI definition");
    let reference_path = map_matches
        .value_of("reference")
        .expect("Presence is ensured by CLI definition");
    let out_file_path = map_matches
        .value_of("output")
        .expect("Presence is ensured by CLI definition");

    let num_threads = value_t!(map_matches.value_of("num_threads"), usize)
        .expect("Presence is ensured by CLI definition");
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Thread pool should not have been initialized before");

    let alignment_parameters = build_alignment_parameters(map_matches);

    if let Err(e) = if map_matches.is_present("dispatcher") {
        info!("Dispatcher mode");
        let port = value_t!(map_matches.value_of("port"), u16).unwrap_or_else(|e| e.exit());
        dispatcher::Dispatcher::new(
            reads_path,
            reference_path,
            out_file_path,
            &alignment_parameters,
        )
        .and_then(|mut dispatcher| dispatcher.run(port))
    } else {
        map::run(
            reads_path,
            reference_path,
            out_file_path,
            &alignment_parameters,
        )
    } {
        error!("{}", e);
    }
}

fn start_worker(arg_matches: &ArgMatches) {
    let host = arg_matches
        .value_of("host")
        .expect("Presence is ensured by CLI definition");
    let port = value_t!(arg_matches.value_of("port"), u16).unwrap_or_else(|e| e.exit());

    let num_threads = value_t!(arg_matches.value_of("num_threads"), usize)
        .expect("Presence is ensured by CLI definition");
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Thread pool should not have been initialized before");

    if let Err(e) = worker::Worker::new(host, port).and_then(|mut worker| worker.run()) {
        error!("{}", e);
    }
}

fn build_alignment_parameters(arg_matches: &ArgMatches) -> AlignmentParameters {
    let library_prep = match arg_matches
        .value_of("library")
        .expect("Presence is ensured by CLI definition")
    {
        "single_stranded" => LibraryPrep::SingleStranded {
            five_prime_overhang: value_t!(arg_matches.value_of("five_prime_overhang"), f32)
                .unwrap_or_else(|e| e.exit()),
            three_prime_overhang: value_t!(arg_matches.value_of("three_prime_overhang"), f32)
                .unwrap_or_else(|e| e.exit()),
        },
        "double_stranded" => LibraryPrep::DoubleStranded(
            value_t!(arg_matches.value_of("five_prime_overhang"), f32).unwrap_or_else(|e| e.exit()),
        ),
        _ => unreachable!(),
    };

    let divergence = value_t!(arg_matches.value_of("divergence"), f32).unwrap_or_else(|e| e.exit());

    let difference_model = SimpleAncientDnaModel::new(
        library_prep,
        value_t!(arg_matches.value_of("ds_deamination_rate"), f32).unwrap_or_else(|e| e.exit()),
        value_t!(arg_matches.value_of("ss_deamination_rate"), f32).unwrap_or_else(|e| e.exit()),
        // Divergence is divided by three because it is used for testing each of the three possible substitutions
        divergence / 3.0,
        arg_matches.is_present("ignore_base_quality"),
    );

    let mismatch_bound = if arg_matches.is_present("poisson_prob") {
        Discrete::new(
            value_t!(arg_matches.value_of("poisson_prob"), f32).unwrap_or_else(|e| e.exit()),
            divergence,
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    } else {
        Continuous::new(
            value_t!(arg_matches.value_of("as_cutoff"), f32).unwrap_or_else(|e| e.exit()) * -1.0,
            value_t!(arg_matches.value_of("as_cutoff_exponent"), f32).unwrap_or_else(|e| e.exit()),
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    };

    AlignmentParameters {
        penalty_gap_open: value_t!(arg_matches.value_of("indel_rate"), f32)
            .unwrap_or_else(|e| e.exit())
            .log2(),
        penalty_gap_extend: value_t!(arg_matches.value_of("gap_extension_penalty"), f32)
            .unwrap_or_else(|e| e.exit())
            * difference_model.get_representative_mismatch_penalty(),
        difference_model: difference_model.into(),
        chunk_size: value_t!(arg_matches.value_of("chunk_size"), usize)
            .unwrap_or_else(|e| e.exit()),
        mismatch_bound,
        gap_dist_ends: value_t!(arg_matches.value_of("gap_dist_ends"), u8)
            .unwrap_or_else(|e| e.exit()),
        stack_limit_abort: arg_matches.is_present("stack_limit_abort"),
    }
}
