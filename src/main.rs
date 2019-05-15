use clap::{
    crate_description, crate_name, crate_version, value_t_or_exit, App, AppSettings, Arg,
    SubCommand,
};

use thrust::sequence_difference_models::{SequenceDifferenceModel, VindijaPWM};
use thrust::{index, map, utils};

fn main() {
    let matches = App::new(crate_name!())
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
                .about("Indexes a genome file")  // TODO
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
                .about("Maps reads to an indexed genome.")  // TODO
                .version(crate_version!())
                .arg(
                    Arg::with_name("reads")
                        .required(true)
                        .short("r")
                        .long("reads")
                        .help("FASTQ file containing adapter-trimmed and quality-controlled reads")
                        .value_name("FASTQ FILE"),
                )
                .arg(
                    Arg::with_name("reference")
                        .required(true)
                        .short("g")
                        .long("reference")
                        .help("Filename of the genome file")
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
                        .short("p")
                        .conflicts_with("max_diff")
                        .default_value("0.04")
                        .help("Minimum probability of the number of mismatches under 0.02 base error rate")
                        .value_name("PROBABILITY")
                        .validator(|v| {
                            let error_message =
                                String::from("Please specify a probability between 0 and 1");
                            let v: f32 = match v.parse() {
                                Ok(s) => s,
                                Err(_) => return Err(error_message),
                            };
                            if (v >= 0.0) && (v <= 1.0) {
                                return Ok(());
                            }
                            Err(error_message)
                        }),
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
            let difference_model = VindijaPWM::new();
            let alignment_parameters = utils::AlignmentParameters {
                base_error_rate: 0.02,
                poisson_threshold: value_t_or_exit!(map_matches.value_of("poisson_prob"), f64),
                penalty_gap_open: 3.5 * difference_model.get_representative_mismatch_penalty(),
                penalty_gap_extend: 1.5 * difference_model.get_representative_mismatch_penalty(),
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
