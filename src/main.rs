use clap::{crate_description, crate_version, Arg, ArgMatches, Command};
use log::{error, info, warn};
#[cfg(all(target_env = "musl"))]
use mimalloc::MiMalloc;
use simple_logger::SimpleLogger;
use time::UtcOffset;

use mapad::{
    distributed::{dispatcher, worker},
    index::indexing,
    map::{
        mapping,
        mismatch_bounds::{Continuous, Discrete},
        sequence_difference_models::{LibraryPrep, SequenceDifferenceModel, SimpleAncientDnaModel},
        AlignmentParameters,
    },
    CRATE_NAME,
};

// Use mimalloc only for musl target
#[cfg(all(target_env = "musl"))]
#[global_allocator]
static ALLOC: MiMalloc = MiMalloc;

fn main() {
    handle_arguments(define_cli());
}

fn define_cli() -> ArgMatches {
    let probability_validator = |v: &str| {
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

    Command::new(CRATE_NAME)
        .about(crate_description!())
        .version(crate_version!())
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("v")
                .global(true)
                .short('v')
                .multiple_occurrences(true)
                .help("Sets the level of verbosity"),
        )
        .arg(
            Arg::new("num_threads")
                .global(true)
                .long("threads")
                .help(&*format!("Maximum number of threads. If 0 or unspecified, {} will select the number of threads automatically.", CRATE_NAME))
                .takes_value(true)
                .default_value("0")
                .value_name("INT")
        )
        .arg(
            Arg::new("port")
                .global(true)
                .long("port")
                .help("TCP port to communicate over")
                .takes_value(true)
                .default_value("3130")
                .value_name("INT")
        )
        .arg(
            Arg::new("seed")
                .global(true)
                .long("seed")
                .help("Seed for the random number generator")
                .takes_value(true)
                .default_value("1234")
                .value_name("INT"),
        )
        .subcommand(
            Command::new("index")
                .about("Indexes a genome file")
                .version(crate_version!())
                .arg(
                    Arg::new("reference")
                        .required(true)
                        .short('g')
                        .long("reference")
                        .help("FASTA file containing the genome to be indexed")
                        .takes_value(true)
                        .value_name("FASTA FILE"),
                )
        )
        .subcommand(
            Command::new("map")
                .about("Maps reads to an indexed genome")
                .version(crate_version!())
                .arg(
                    Arg::new("reads")
                        .required(true)
                        .short('r')
                        .long("reads")
                        .help("BAM or FASTQ file that contains the reads to be aligned")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("reference")
                        .required(true)
                        .short('g')
                        .long("reference")
                        .help("Prefix of the file names of the index files. The reference FASTA file itself does not need to be present.")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output")
                        .required(true)
                        .short('o')
                        .long("output")
                        .help("Path to output BAM file")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("poisson_prob")
                        .short('p')
                        .group("allowed_mm")
                        .help("Minimum probability of the number of mismatches under `-D` base error rate")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("as_cutoff")
                        .short('c')
                        .group("allowed_mm")
                        .help("Per-base average alignment score cutoff (-c > AS / read_len^e ?)")
                        .takes_value(true)
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("as_cutoff_exponent")
                        .short('e')
                        .help("Exponent to be applied to the read length (ignored if `-c` is not used)")
                        .takes_value(true)
                        .default_value("1.0")
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("library")
                        .required(true)
                        .short('l')
                        .long("library")
                        .help("Library preparation method")
                        .takes_value(true)
                        .possible_values(&["single_stranded", "double_stranded"])
                        .value_name("STRING")
                )
                .arg(
                    Arg::new("five_prime_overhang")
                        .required(true)
                        .short('f')
                        .help("5'-overhang length parameter")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("three_prime_overhang")
                        .required_if_eq("library", "single_stranded")
                        .short('t')
                        .help("3'-overhang length parameter")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("ds_deamination_rate")
                        .required(true)
                        .short('d')
                        .help("Deamination rate in double-stranded stem of a read")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("ss_deamination_rate")
                        .required(true)
                        .short('s')
                        .help("Deamination rate in single-stranded ends of a read")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("divergence")
                        .short('D')
                        .help("Divergence / base error rate")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .default_value("0.02")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("indel_rate")
                        .required(true)
                        .short('i')
                        .help("Expected rate of indels between reads and reference")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("gap_extension_penalty")
                        .short('x')
                        .help("Gap extension penalty as a fraction of the representative mismatch penalty")
                        .takes_value(true)
                        .value_name("FLOAT")
                        .default_value("1.0")
                        .validator(probability_validator),
                )
                .arg(
                    Arg::new("chunk_size")
                        .long("batch_size")
                        .help("The number of reads that are processed in parallel")
                        .takes_value(true)
                        .default_value("250000")
                        .value_name("INT")
                )
                .arg(
                    Arg::new("ignore_base_quality")
                        .long("ignore_base_quality")
                        .help("Ignore base qualities in scoring models")
                        .takes_value(false)
                )
                .arg(
                    Arg::new("dispatcher")
                        .long("dispatcher")
                        .help("Run in dispatcher mode for distributed computing in a network. Needs workers to be spawned externally to distribute work among them.")
                        .takes_value(false)
                )
                .arg(
                    Arg::new("gap_dist_ends")
                        .long("gap_dist_ends")
                        .help("Disallow gaps at read ends (configurable range)")
                        .takes_value(true)
                        .default_value("5")
                        .value_name("INT")
                )
                .arg(
                    Arg::new("stack_limit_abort")
                        .long("stack_limit_abort")
                        .help("Abort alignment when stack size limit is reached instead of trying to recover.")
                        .takes_value(false)
                ),

        )
        .subcommand(
            Command::new("worker")
                .about("Spawns worker")
                .version(crate_version!())
                .arg(
                    Arg::new("host")
                        .required(true)
                        .long("host")
                        .help("Hostname or IP address of the running dispatcher node")
                        .takes_value(true)
                        .value_name("HOST")
                )
        )
        .get_matches()
}

fn handle_arguments(matches: ArgMatches) {
    let mut logger_builder = SimpleLogger::new().with_level(match matches.occurrences_of("v") {
        0 => log::LevelFilter::Info,
        1 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    });
    if let Ok(local_offset) = UtcOffset::current_local_offset() {
        logger_builder = logger_builder.with_utc_offset(local_offset);
    } else {
        warn!("Use UTC timestamps for log entries");
        logger_builder = logger_builder.with_utc_timestamps();
    }
    logger_builder.init().expect("This is not expected to fail");

    let seed = matches.value_of_t("seed").unwrap_or_else(|e| e.exit());
    match matches.subcommand() {
        Some(("index", arg_matches)) => {
            start_indexer(arg_matches, seed);
        }
        Some(("map", arg_matches)) => {
            start_mapper(arg_matches, seed);
        }
        Some(("worker", arg_matches)) => {
            start_worker(arg_matches, seed);
        }
        _ => unreachable!(),
    }
}

fn start_indexer(arg_matches: &ArgMatches, seed: u64) {
    let reference_path = arg_matches
        .value_of("reference")
        .expect("Presence is ensured by CLI definition");

    if let Err(e) = indexing::run(reference_path, seed) {
        error!("{}", e);
    }
}

fn start_mapper(map_matches: &ArgMatches, _seed: u64) {
    let reads_path = map_matches
        .value_of("reads")
        .expect("Presence is ensured by CLI definition");
    let reference_path = map_matches
        .value_of("reference")
        .expect("Presence is ensured by CLI definition");
    let out_file_path = map_matches
        .value_of("output")
        .expect("Presence is ensured by CLI definition");

    let num_threads = map_matches
        .value_of_t("num_threads")
        .unwrap_or_else(|e| e.exit());
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Thread pool should not have been initialized before");

    let alignment_parameters = build_alignment_parameters(map_matches);

    if let Err(e) = if map_matches.is_present("dispatcher") {
        info!("Dispatcher mode");
        let port = map_matches.value_of_t("port").unwrap_or_else(|e| e.exit());
        dispatcher::Dispatcher::new(
            reads_path,
            reference_path,
            out_file_path,
            &alignment_parameters,
        )
        .and_then(|mut dispatcher| dispatcher.run(port))
    } else {
        mapping::run(
            reads_path,
            reference_path,
            out_file_path,
            &alignment_parameters,
        )
    } {
        error!("{}", e);
    }
}

fn start_worker(arg_matches: &ArgMatches, _seed: u64) {
    let host = arg_matches
        .value_of("host")
        .expect("Presence is ensured by CLI definition");
    let port = arg_matches.value_of_t("port").unwrap_or_else(|e| e.exit());

    let num_threads = arg_matches
        .value_of_t("num_threads")
        .unwrap_or_else(|e| e.exit());
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
            five_prime_overhang: arg_matches
                .value_of_t("five_prime_overhang")
                .unwrap_or_else(|e| e.exit()),
            three_prime_overhang: arg_matches
                .value_of_t("three_prime_overhang")
                .unwrap_or_else(|e| e.exit()),
        },
        "double_stranded" => LibraryPrep::DoubleStranded(
            arg_matches
                .value_of_t("five_prime_overhang")
                .unwrap_or_else(|e| e.exit()),
        ),
        _ => unreachable!(),
    };

    let divergence = arg_matches
        .value_of_t("divergence")
        .unwrap_or_else(|e| e.exit());

    let difference_model = SimpleAncientDnaModel::new(
        library_prep,
        arg_matches
            .value_of_t("ds_deamination_rate")
            .unwrap_or_else(|e| e.exit()),
        arg_matches
            .value_of_t("ss_deamination_rate")
            .unwrap_or_else(|e| e.exit()),
        // Divergence is divided by three because it is used for testing each of the three possible substitutions
        divergence / 3.0,
        arg_matches.is_present("ignore_base_quality"),
    );

    let mismatch_bound = if arg_matches.is_present("poisson_prob") {
        Discrete::new(
            arg_matches
                .value_of_t("poisson_prob")
                .unwrap_or_else(|e| e.exit()),
            divergence,
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    } else {
        Continuous::new(
            arg_matches
                .value_of_t::<f32>("as_cutoff")
                .unwrap_or_else(|e| e.exit())
                * -1.0,
            arg_matches
                .value_of_t("as_cutoff_exponent")
                .unwrap_or_else(|e| e.exit()),
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    };

    AlignmentParameters {
        penalty_gap_open: arg_matches
            .value_of_t::<f32>("indel_rate")
            .unwrap_or_else(|e| e.exit())
            .log2(),
        penalty_gap_extend: arg_matches
            .value_of_t::<f32>("gap_extension_penalty")
            .unwrap_or_else(|e| e.exit())
            * difference_model.get_representative_mismatch_penalty(),
        difference_model: difference_model.into(),
        chunk_size: arg_matches
            .value_of_t("chunk_size")
            .unwrap_or_else(|e| e.exit()),
        mismatch_bound,
        gap_dist_ends: arg_matches
            .value_of_t("gap_dist_ends")
            .unwrap_or_else(|e| e.exit()),
        stack_limit_abort: arg_matches.is_present("stack_limit_abort"),
    }
}
