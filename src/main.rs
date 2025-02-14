use clap::{crate_authors, crate_description, value_parser, Arg, ArgAction, ArgMatches, Command};
use log::{error, info, warn};
use simple_logger::SimpleLogger;
use time::UtcOffset;

use mapad::{
    build_info,
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
#[cfg(target_env = "musl")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    handle_arguments(&define_cli());
}

fn define_cli() -> ArgMatches {
    /// Parse `f32` float and check if the value is between 0 and 1.
    /// `clap` impls `TypedValueParser` for functions with the following signature,
    /// so it is directly pluggable into its `Arg::value_parser()` method.
    fn parse_validate_prob(raw_arg: &str) -> Result<f32, clap::Error> {
        raw_arg
            .parse()
            .ok()
            .and_then(|value| (0.0..=1.0).contains(&value).then_some(value))
            .ok_or_else(|| clap::Error::new(clap::error::ErrorKind::ValueValidation))
    }

    Command::new(CRATE_NAME)
        .about(crate_description!())
        .version(build_info::get_software_version())
        .author(crate_authors!())
        .subcommand_required(true)
        .arg_required_else_help(true)
        .arg(
            Arg::new("v")
                .global(true)
                .short('v')
                .action(ArgAction::Count)
                .help("Sets the level of verbosity")
        )
        .arg(
            Arg::new("num_threads")
                .global(true)
                .long("threads")
                .help(format!("Maximum number of threads. If 0, {CRATE_NAME} will select the number of threads automatically."))
                .default_value("1")
                .value_name("INT")
                .value_parser(value_parser!(usize))
        )
        .arg(
            Arg::new("port")
                .global(true)
                .long("port")
                .help("TCP port to communicate over")
                .default_value("3130")
                .value_name("INT")
                .value_parser(value_parser!(u16))
        )
        .arg(
            Arg::new("seed")
                .global(true)
                .long("seed")
                .help("Seed for the random number generator")
                .default_value("1234")
                .value_name("INT")
                .value_parser(value_parser!(u64))
        )
        .subcommand(
            Command::new("index")
                .about("Indexes a genome file")
                .arg(
                    Arg::new("reference")
                        .required(true)
                        .short('g')
                        .long("reference")
                        .help("FASTA file containing the genome to be indexed")
                        .value_name("FASTA FILE")
                )
        )
        .subcommand(
            Command::new("map")
                .about("Maps reads to an indexed genome")
                .arg(
                    Arg::new("reads")
                        .required(true)
                        .short('r')
                        .long("reads")
                        .help("BAM/CRAM/FASTQ or FASTQ.GZ file that contains the reads to be aligned. Specify \"-\" for reading from stdin.")
                        .value_name("STRING")
                )
                .arg(
                    Arg::new("reference")
                        .required(true)
                        .short('g')
                        .long("reference")
                        .help("Prefix of the file names of the index files. The reference FASTA file itself does not need to be present.")
                        .value_name("STRING")
                )
                .arg(
                    Arg::new("output")
                        .required(true)
                        .short('o')
                        .long("output")
                        .help("Path to output BAM file")
                        .value_name("STRING")
                )
                .arg(
                    Arg::new("poisson_prob")
                        .short('p')
                        .help("Minimum probability of the number of mismatches under `-D` base error rate")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                        .conflicts_with("as_cutoff")
                        .required_unless_present("as_cutoff")
                )
                .arg(
                    Arg::new("as_cutoff")
                        .short('c')
                        .help("Per-base average alignment score cutoff (-c > AS / read_len^e ?)")
                        .value_name("FLOAT")
                        .value_parser(value_parser!(f32))
                )
                .arg(
                    Arg::new("as_cutoff_exponent")
                        .short('e')
                        .help("Exponent to be applied to the read length (ignored if `-c` is not used)")
                        .default_value("1.0")
                        .value_name("FLOAT")
                    .value_parser(value_parser!(f32))
                )
                .arg(
                    Arg::new("library")
                        .required(true)
                        .short('l')
                        .long("library")
                        .help("Library preparation method")
                        .value_parser(["single_stranded", "double_stranded"])
                        .value_name("STRING")
                )
                .arg(
                    Arg::new("five_prime_overhang")
                        .required(true)
                        .short('f')
                        .help("5'-overhang length parameter")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("three_prime_overhang")
                        .required_if_eq("library", "single_stranded")
                        .short('t')
                        .help("3'-overhang length parameter")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("ds_deamination_rate")
                        .required(true)
                        .short('d')
                        .help("Deamination rate in double-stranded stem of a read")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("ss_deamination_rate")
                        .required(true)
                        .short('s')
                        .help("Deamination rate in single-stranded ends of a read")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("divergence")
                        .short('D')
                        .help("Divergence / base error rate")
                        .value_name("FLOAT")
                        .default_value("0.02")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("indel_rate")
                        .required(true)
                        .short('i')
                        .help("Expected rate of indels between reads and reference")
                        .value_name("FLOAT")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("gap_extension_penalty")
                        .short('x')
                        .help("Gap extension penalty as a fraction of the representative mismatch penalty")
                        .value_name("FLOAT")
                        .default_value("1.0")
                        .value_parser(parse_validate_prob)
                )
                .arg(
                    Arg::new("chunk_size")
                        .long("batch_size")
                        .help("The number of reads that are processed in parallel")
                        .default_value("250000")
                        .value_name("INT")
                        .value_parser(value_parser!(usize))
                )
                .arg(
                    Arg::new("ignore_base_quality")
                        .long("ignore_base_quality")
                        .help("Ignore base qualities in scoring models")
                        .action(ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("dispatcher")
                        .long("dispatcher")
                        .help("Run in dispatcher mode for distributed computing in a network. Needs workers to be spawned externally to distribute work among them.")
                        .action(ArgAction::SetTrue)
                )
                .arg(
                    Arg::new("gap_dist_ends")
                        .long("gap_dist_ends")
                        .help("Disallow gaps at read ends (configurable range)")
                        .default_value("5")
                        .value_name("INT")
                        .value_parser(value_parser!(u8))
                )
                .arg(
                    Arg::new("max_num_gaps_open")
                        .long("max_num_gaps_open")
                        .help("Max. number of opened gaps")
                        .default_value("2")
                        .value_name("INT")
                        .value_parser(value_parser!(u8))
                )
                .arg(
                    Arg::new("search_limit_recovery")
                        .long("search_limit_recovery")
                        .help(format!(
                            "Mapping of reads which are particularly difficult to align can exceed \
                            internal search space limits. With this option enabled, {CRATE_NAME} \
                            will try to recover from these cases by discarding low-scoring \
                            sub-alignments instead of reporting the read as unmapped. Enabling \
                            this option will potentially slow down the mapping."
                        ))
                        .action(ArgAction::SetFalse),
                )
                .arg(
                    Arg::new("force_overwrite")
                        .long("force_overwrite")
                        .help(format!("Force {CRATE_NAME} to overwrite the output BAM file."))
                        .action(ArgAction::SetTrue)
                )
        )
        .subcommand(
            Command::new("worker")
                .about("Spawns worker")
                .arg(
                    Arg::new("host")
                        .required(true)
                        .long("host")
                        .help("Hostname or IP address of the running dispatcher node")
                        .value_name("HOST")
                )
        )
        .get_matches()
}

fn handle_arguments(matches: &ArgMatches) {
    let mut logger_builder = SimpleLogger::new().with_level(match matches.get_count("v") {
        0 => log::LevelFilter::Info,
        1 => log::LevelFilter::Debug,
        2.. => log::LevelFilter::Trace,
    });
    if let Ok(local_offset) = UtcOffset::current_local_offset() {
        logger_builder = logger_builder.with_utc_offset(local_offset);
    } else {
        warn!("Use UTC timestamps for log entries");
        logger_builder = logger_builder.with_utc_timestamps();
    }
    logger_builder.init().expect("This is not expected to fail");

    let seed = *matches
        .get_one("seed")
        .expect("Presence ensured by CLI definition");
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
        .get_one::<String>("reference")
        .expect("Presence is ensured by CLI definition");

    if let Err(e) = indexing::run(reference_path, seed) {
        error!("{}", e);
    }
}

fn start_mapper(map_matches: &ArgMatches, _seed: u64) {
    let reads_path = map_matches
        .get_one::<String>("reads")
        .expect("Presence is ensured by CLI definition");
    let reference_path = map_matches
        .get_one::<String>("reference")
        .expect("Presence is ensured by CLI definition");
    let out_file_path = map_matches
        .get_one::<String>("output")
        .expect("Presence is ensured by CLI definition");
    let force_overwrite = map_matches.get_flag("force_overwrite");

    let num_threads = *map_matches
        .get_one("num_threads")
        .expect("Presence ensured by CLI definition");
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Thread pool should not have been initialized before");

    let alignment_parameters = build_alignment_parameters(map_matches);

    if let Err(e) = if map_matches.get_flag("dispatcher") {
        info!("Dispatcher mode");
        let port = *map_matches
            .get_one("port")
            .expect("Presence ensured by CLI definition");
        dispatcher::Dispatcher::new(
            reads_path,
            reference_path,
            out_file_path,
            force_overwrite,
            &alignment_parameters,
        )
        .and_then(|mut dispatcher| dispatcher.run(port))
    } else {
        mapping::run(
            reads_path,
            reference_path,
            out_file_path,
            force_overwrite,
            &alignment_parameters,
        )
    } {
        error!("{}", e);
    }
}

fn start_worker(arg_matches: &ArgMatches, _seed: u64) {
    let host = arg_matches
        .get_one::<String>("host")
        .expect("Presence is ensured by CLI definition");
    let port = *arg_matches
        .get_one("port")
        .expect("Presence ensured by CLI definition");

    let num_threads = *arg_matches
        .get_one("num_threads")
        .expect("Presence ensured by CLI definition");
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Thread pool should not have been initialized before");

    if let Err(e) = worker::Worker::new(host, port).and_then(|mut worker| worker.run()) {
        error!("{}", e);
    }
}

fn build_alignment_parameters(arg_matches: &ArgMatches) -> AlignmentParameters {
    let library_prep = match &**arg_matches
        .get_one::<String>("library")
        .expect("Presence is ensured by CLI definition")
    {
        "single_stranded" => LibraryPrep::SingleStranded {
            five_prime_overhang: *arg_matches
                .get_one("five_prime_overhang")
                .expect("Presence ensured by CLI definition"),
            three_prime_overhang: *arg_matches
                .get_one("three_prime_overhang")
                .expect("Presence ensured by CLI definition"),
        },
        "double_stranded" => LibraryPrep::DoubleStranded(
            *arg_matches
                .get_one("five_prime_overhang")
                .expect("Presence ensured by CLI definition"),
        ),
        _ => unreachable!(),
    };

    let divergence = *arg_matches
        .get_one("divergence")
        .expect("Presence ensured by CLI definition");

    let difference_model = SimpleAncientDnaModel::new(
        library_prep,
        *arg_matches
            .get_one("ds_deamination_rate")
            .expect("Presence ensured by CLI definition"),
        *arg_matches
            .get_one("ss_deamination_rate")
            .expect("Presence ensured by CLI definition"),
        // Divergence is divided by three because it is used for testing each of the three possible substitutions
        divergence / 3.0,
        arg_matches.get_flag("ignore_base_quality"),
    );

    let mismatch_bound = if let Some(&pprob) = arg_matches.get_one("poisson_prob") {
        Discrete::new(
            pprob,
            divergence,
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    } else {
        Continuous::new(
            arg_matches
                .get_one::<f32>("as_cutoff")
                .expect("Presence ensured by CLI definition")
                * -1.0,
            *arg_matches
                .get_one::<f32>("as_cutoff_exponent")
                .expect("Presence ensured by CLI definition"),
            difference_model.get_representative_mismatch_penalty(),
        )
        .into()
    };

    AlignmentParameters {
        penalty_gap_open: arg_matches
            .get_one::<f32>("indel_rate")
            .expect("Presence ensured by CLI definition")
            .log2(),
        penalty_gap_extend: arg_matches
            .get_one::<f32>("gap_extension_penalty")
            .expect("Presence ensured by CLI definition")
            * difference_model.get_representative_mismatch_penalty(),
        difference_model: difference_model.into(),
        chunk_size: *arg_matches
            .get_one("chunk_size")
            .expect("Presence ensured by CLI definition"),
        mismatch_bound,
        gap_dist_ends: *arg_matches
            .get_one("gap_dist_ends")
            .expect("Presence ensured by CLI definition"),
        stack_limit_abort: arg_matches.get_flag("stack_limit_abort"),
        max_num_gaps_open: *arg_matches
            .get_one("max_num_gaps_open")
            .expect("Presence ensured by CLI definition"),
    }
}
