use std::{
    borrow::Cow,
    cell::RefCell,
    collections::BinaryHeap,
    env,
    fs::OpenOptions,
    io::{self, Write},
    num::NonZeroUsize,
    path::Path,
    str,
    time::{Duration, Instant},
};

use bio::{alphabets::dna, data_structures::suffix_array::SuffixArray};
use bstr::{BString, ByteSlice};
use clap::crate_description;
use log::{debug, info, trace};
use min_max_heap::MinMaxHeap;
use noodles::{
    bam,
    sam::{self, alignment::io::Write as SamWrite},
};
use rand::RngCore;
use rayon::prelude::*;

use crate::{
    build_info,
    errors::{Error, Result},
    index::{
        load_id_pos_map_from_path, load_index_from_path, load_original_symbols_from_path,
        load_suffix_array_from_path, FastaIdPositions, OriginalSymbols,
    },
    map::{
        backtrack_tree::Tree,
        bi_d_array::BiDArray,
        fmd_index::RtFmdIndex,
        input_chunk_reader::{InputSource, TaskQueue},
        mismatch_bounds::{MismatchBound, MismatchBoundDispatch},
        prrange::PrRange,
        record::{extract_edit_operations, EditOperation, Record},
        sequence_difference_models::{SequenceDifferenceModel, SequenceDifferenceModelDispatch},
        AlignmentParameters, AlternativeAlignments, Direction, GapState, HitInterval,
        IntToCoordOutput, MismatchSearchStackFrame,
    },
    CRATE_NAME,
};

// These settings lead to a memory consumption of ~245 MiB per thread
pub const STACK_LIMIT: u32 = 2_000_000;
pub const EDIT_TREE_LIMIT: u32 = 10_000_000;

/// Loads index files and launches the mapping process
pub fn run(
    reads_path: &str,
    reference_path: &str,
    out_file_path: &str,
    force_overwrite: bool,
    alignment_parameters: &AlignmentParameters,
) -> Result<()> {
    let reads_path = Path::new(reads_path);
    let out_file_path = Path::new(out_file_path);
    if !reads_path.exists() && reads_path.to_str() != Some("-") {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "The given input file could not be found",
        )
        .into());
    }

    info!("Load FMD-index");
    let fmd_index = load_index_from_path(reference_path)?;

    info!("Load suffix array");
    let sampled_suffix_array_owned = load_suffix_array_from_path(reference_path)?;
    let suffix_array = sampled_suffix_array_owned.into_sampled_suffix_array(
        &fmd_index.bwt,
        &fmd_index.less,
        &fmd_index.occ,
    );

    info!("Load position map");
    let identifier_position_map = load_id_pos_map_from_path(reference_path)?;

    info!("Load original symbols");
    let original_symbols = load_original_symbols_from_path(reference_path)?;

    let mut out_file = bam::io::Writer::new(
        OpenOptions::new()
            .read(false)
            .write(true)
            // If .create_new(true) is set, .create() and .truncate() are ignored
            .create_new(!force_overwrite)
            .create(true)
            .truncate(true)
            .open(out_file_path)?,
    );

    info!("Map reads");
    let mut input_source = InputSource::from_path(reads_path)?;
    let out_header = create_bam_header(input_source.header(), &identifier_position_map)?;
    out_file.write_header(&out_header)?;
    run_inner(
        input_source.task_queue(alignment_parameters.chunk_size),
        &fmd_index,
        &suffix_array,
        alignment_parameters,
        &identifier_position_map,
        &original_symbols,
        &out_header,
        &mut out_file,
    )?;

    info!("Done");
    Ok(())
}

/// This part has been extracted from the main run() function to allow static dispatch based on the
/// input file type
#[allow(clippy::too_many_arguments)]
fn run_inner<S, T, W>(
    records: TaskQueue<T>,
    fmd_index: &RtFmdIndex,
    suffix_array: &S,
    alignment_parameters: &AlignmentParameters,
    identifier_position_map: &FastaIdPositions,
    original_symbols: &OriginalSymbols,
    out_header: &sam::Header,
    out_file: &mut bam::io::Writer<W>,
) -> Result<()>
where
    S: SuffixArray + Send + Sync,
    T: Iterator<Item = Result<Record>>,
    W: Write,
{
    thread_local! {
        static STACK_BUF: RefCell<MinMaxHeap<MismatchSearchStackFrame>> = RefCell::new(MinMaxHeap::with_capacity(STACK_LIMIT as usize + 9));
        static TREE_BUF: RefCell<Tree<EditOperation>> = RefCell::new(Tree::with_capacity(EDIT_TREE_LIMIT + 9));
    }

    for chunk in records {
        debug!("Map chunk of records");
        let bam_records = chunk
            .get_records()
            .into_par_iter()
            .map(|record| {
                STACK_BUF.with(|stack_buf| {
                    TREE_BUF.with(|tree_buf| {
                        let start = Instant::now();
                        // Here we call the instances of the generic `k_mismatch_search` function. Currently, only
                        // `SimpleAncientDnaModel` is used, the others merely serve as an example.
                        let hit_intervals = match &alignment_parameters.difference_model {
                            SequenceDifferenceModelDispatch::SimpleAncientDnaModel(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                            SequenceDifferenceModelDispatch::TestDifferenceModel(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                            SequenceDifferenceModelDispatch::VindijaPwm(sdm) => {
                                match &alignment_parameters.mismatch_bound {
                                    MismatchBoundDispatch::Discrete(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::Continuous(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                    MismatchBoundDispatch::TestBound(mb) => k_mismatch_search(
                                        &record.sequence,
                                        &record.base_qualities,
                                        alignment_parameters,
                                        fmd_index,
                                        &mut stack_buf.borrow_mut(),
                                        &mut tree_buf.borrow_mut(),
                                        sdm,
                                        mb,
                                    ),
                                }
                            }
                        };
                        let duration = start.elapsed();

                        (record, hit_intervals, duration)
                    })
                })
            })
            .map_init(
                rand::thread_rng,
                |mut rng, (record, hit_interval, duration)| -> Result<sam::alignment::RecordBuf> {
                    intervals_to_bam(
                        record,
                        hit_interval,
                        suffix_array,
                        identifier_position_map,
                        original_symbols,
                        Some(&duration),
                        alignment_parameters,
                        &mut rng,
                    )
                },
            )
            .collect::<Result<Vec<_>>>()?;

        debug!("Write chunk of BAM records to output file");
        for record in &bam_records {
            out_file.write_alignment_record(out_header, record)?;
        }
    }
    Ok(())
}

/// Creates basic BAM header with mandatory fields pre-populated.
/// If provided, it will copy most data from an input BAM header (`@SQ` lines are removed).
pub fn create_bam_header(
    src_header: Option<&sam::Header>,
    identifier_position_map: &FastaIdPositions,
) -> Result<sam::Header> {
    use sam::header::record::value::map::{self, Map};

    let mut sam_header_builder = sam::Header::builder();
    let mut sam_header_header = Map::<map::Header>::new(map::header::Version::new(1, 6));
    sam_header_header.other_fields_mut().insert(
        map::header::tag::SORT_ORDER,
        map::header::sort_order::UNSORTED.into(),
    );

    let mut program_id = Cow::from(CRATE_NAME);
    let program_builder = {
        let cmdline = {
            let mut out = env::args().fold(String::new(), |acc, part| acc + &part + " ");
            let _ = out.pop();
            out
        };
        Map::<map::Program>::builder()
            .insert(map::program::tag::NAME, CRATE_NAME)
            .insert(
                map::program::tag::VERSION,
                build_info::get_software_version(),
            )
            .insert(map::program::tag::DESCRIPTION, crate_description!())
            .insert(map::program::tag::COMMAND_LINE, cmdline)
    };

    if let Some(src_header) = src_header {
        // We've got a header to copy some data from

        // @PG chain of old entries
        for (id, pg) in src_header.programs().as_ref() {
            sam_header_builder = sam_header_builder.add_program(id.as_bstr(), pg.to_owned());
        }

        // Ensure @PG/ID is unique
        {
            let pg_id_count = src_header
                .programs()
                .as_ref()
                .keys()
                .filter(|id| {
                    id.as_bstr() == program_id
                        || id.starts_with(format!("{program_id}.").as_bytes())
                })
                .count();
            if pg_id_count > 0 {
                program_id = Cow::from(format!("{program_id}.{pg_id_count}"));
            }
        }

        for comment in src_header.comments() {
            sam_header_builder = sam_header_builder.add_comment(comment.as_bstr());
        }

        for (id, read_group) in src_header.read_groups() {
            sam_header_builder =
                sam_header_builder.add_read_group(id.as_bstr(), read_group.clone());
        }
    }

    // @SQ entries
    for identifier_position in identifier_position_map.iter() {
        sam_header_builder = sam_header_builder.add_reference_sequence(
            identifier_position.identifier.clone(),
            Map::<map::ReferenceSequence>::new(
                usize::try_from(identifier_position.end - identifier_position.start + 1).and_then(NonZeroUsize::try_from).map_err(
                    |_e| {
                        Error::InvalidIndex(
                            "Could not create header. Reference genome size seems too large to be processed on this machine.".into(),
                        )
                    },
                )?,
            )
        );
    }

    // Finalize header
    sam_header_builder = sam_header_builder.set_header(sam_header_header);
    let mut sam_header = sam_header_builder.build();

    // Finalize and add new @PG entry
    let program = program_builder
        .build()
        .expect("@PG entry is generated internally");
    sam_header
        .programs_mut()
        .add(program_id.as_ref(), program)?;

    Ok(sam_header)
}

/// Convert suffix array intervals to positions and BAM records
#[allow(clippy::too_many_arguments)]
pub fn intervals_to_bam<R, S>(
    input_record: Record,
    intervals: BinaryHeap<HitInterval>,
    suffix_array: &S,
    identifier_position_map: &FastaIdPositions,
    original_symbols: &OriginalSymbols,
    duration: Option<&Duration>,
    alignment_parameters: &AlignmentParameters,
    rng: &mut R,
) -> Result<sam::alignment::RecordBuf>
where
    R: RngCore,
    S: SuffixArray,
{
    let hits_found = !intervals.is_empty();

    let mut intervals = intervals.into_sorted_vec();

    while let Some(best_alignment) = intervals.pop() {
        // REMEMBER: The best-scoring alignment has been `pop()`ed,
        // so `intervals`, from here, only contains suboptimal alignments
        let mut best_i2co_iter =
            interval2coordinate(&best_alignment, suffix_array, identifier_position_map, rng)?;

        // The first position is reported in the BAM record
        match best_i2co_iter.next() {
            Some(best_i2co) => {
                let updated_best_alignment_interval_size =
                    best_alignment.interval.size - best_i2co.num_skipped;

                // Alternative hits (same or worse score) are reported via BAM auxiliary fields
                let alternative_hits = {
                    // Chain best multi-mapping hits and suboptimal alignments
                    let xa = best_i2co_iter
                        .chain(
                            intervals
                                .iter()
                                .rev()
                                // Remove duplicates (coming from e.g. InDels in homopolymers)
                                .filter(|suboptimal_alignment| {
                                    suboptimal_alignment.interval != best_alignment.interval
                                })
                                // Report all valid positions of each suboptimal alignment
                                .filter_map(|suboptimal_alignment| {
                                    interval2coordinate(
                                        suboptimal_alignment,
                                        suffix_array,
                                        identifier_position_map,
                                        rng,
                                    )
                                    .ok()
                                })
                                .flatten(),
                        )
                        .map(|i2co| {
                            let (pre_cigar, md, nm) = i2co.interval.edit_operations.to_bam_fields(
                                i2co.strand,
                                i2co.absolute_pos,
                                original_symbols,
                            );
                            let cigar_string = {
                                let cigar = pre_cigar
                                    .into_iter()
                                    .collect::<sam::alignment::record_buf::Cigar>();
                                let mut buf = Vec::new();
                                sam::io::writer::record::write_cigar(&mut buf, &cigar)
                                    .expect("Generated internally. Should be correct.");
                                buf
                            };
                            format!(
                                "{},{}{},{},{},{},{},{:.2};",
                                i2co.contig_name,
                                match i2co.strand {
                                    Direction::Forward => '+',
                                    Direction::Backward => '-',
                                },
                                i2co.relative_pos + 1,
                                String::from_utf8_lossy(&cigar_string),
                                String::from_utf8_lossy(&md),
                                nm,
                                i2co.interval.interval.size,
                                i2co.interval.alignment_score,
                            )
                        })
                        .take(2)
                        .collect::<String>();

                    AlternativeAlignments {
                        x0: updated_best_alignment_interval_size
                            .try_into()
                            .unwrap_or(i32::MAX),
                        x1: intervals
                            .iter()
                            .filter(|suboptimal_alignment| {
                                suboptimal_alignment.interval != best_alignment.interval
                            })
                            .map(|suboptimal_alignment| suboptimal_alignment.interval.size)
                            .sum::<usize>()
                            .try_into()
                            .unwrap_or(i32::MAX),
                        xa,
                        xs: intervals
                            .last()
                            .map(|suboptimal_alignment| suboptimal_alignment.alignment_score)
                            .unwrap_or_default(),
                        xt: match updated_best_alignment_interval_size {
                            0 => 'N',
                            1 => 'U',
                            _ => 'R',
                        },
                    }
                };

                return create_bam_record(
                    input_record,
                    Some(best_i2co.relative_pos),
                    Some(best_i2co.absolute_pos),
                    Some(&best_alignment),
                    Some(estimate_mapping_quality(
                        &best_alignment,
                        updated_best_alignment_interval_size,
                        &intervals,
                        alignment_parameters,
                    )),
                    Some(best_i2co.tid),
                    Some(best_i2co.strand),
                    duration,
                    Some(alternative_hits),
                    original_symbols,
                );
            }
            None => {
                debug!("Mapped coordinate overlaps contig boundary. Try again with next best hit.");
            }
        }
    }

    if hits_found {
        debug!(
            "Hits could not be mapped to valid coordinates. Report read \"{}\" as unmapped.",
            input_record
        );
    }

    // No match found, report unmapped read
    create_bam_record(
        input_record,
        None,
        None,
        None,
        Some(0),
        None,
        None,
        duration,
        None,
        original_symbols,
    )
}

/// Computes optimal per-base alignment scores for a read,
/// conditioned on base qualities and scoring model.
/// This function panics if `pattern` and `base_qualities` are not of equal length.
pub fn compute_optimal_scores<SDM>(
    pattern: &[u8],
    base_qualities: &[u8],
    difference_model: &SDM,
) -> Vec<f32>
where
    SDM: SequenceDifferenceModel,
{
    pattern
        .iter()
        .zip(base_qualities)
        .enumerate()
        .map(|(i, (&base, &quality))| {
            difference_model.get_min_penalty(i, pattern.len(), base, quality, false)
        })
        .collect::<Vec<_>>()
}

fn interval2coordinate<'a, S, R>(
    interval: &'a HitInterval,
    suffix_array: &'a S,
    identifier_position_map: &'a FastaIdPositions,
    rng: &mut R,
) -> Result<impl Iterator<Item = IntToCoordOutput<'a>> + 'a>
where
    R: RngCore,
    S: SuffixArray,
{
    // Calculate length of reference strand including sentinel characters
    let strand_len = suffix_array.len() / 2;
    // Get amount of positions in the genome covered by the read
    let effective_read_len = interval.edit_operations.effective_len();

    Ok(usize::try_from(rng.next_u32())
        .ok()
        .and_then(|seed| PrRange::try_from_range(&interval.interval.range_fwd(), seed))
        .ok_or_else(|| {
            Error::InvalidIndex("Could not enumerate possible reference positions".into())
        })?
        .enumerate() // Count the number of skipped invalid positions
        .filter_map(move |(i, sar_pos)| {
            suffix_array
                .get(sar_pos)
                // Determine strand
                .map(|absolute_pos| {
                    if absolute_pos < strand_len {
                        (absolute_pos, Direction::Forward)
                    } else {
                        (
                            // Due to the sentinel of the forward strand, we'd need to subtract 2
                            // to map a reverse strand position to the forward strand.
                            // The following mapping is a slightly shorter version of the equivalent
                            // `suffix_array.len() - absolute_pos - (effective_read_len - 1) - 2`.
                            suffix_array.len() - absolute_pos - effective_read_len - 1,
                            Direction::Backward,
                        )
                    }
                })
                // Convert to relative coordinates
                .and_then(|(absolute_pos, strand)| {
                    if let Some((tid, rel_pos, contig_name)) = identifier_position_map
                        .get_reference_identifier(absolute_pos, effective_read_len)
                    {
                        Some(IntToCoordOutput {
                            tid,
                            interval,
                            contig_name,
                            relative_pos: rel_pos,
                            absolute_pos,
                            strand,
                            num_skipped: i,
                        })
                    } else {
                        None
                    }
                })
        }))
}

/// Estimate mapping quality based on the number of hits for a particular read, its alignment score,
/// and its base qualities
#[allow(clippy::assertions_on_constants)]
fn estimate_mapping_quality(
    best_alignment: &HitInterval,
    best_alignment_interval_size: usize,
    other_alignments: &[HitInterval],
    alignment_parameters: &AlignmentParameters,
) -> u8 {
    const MAX_MAPQ: u8 = 37;
    const MIN_MAPQ_UNIQ: u8 = 20;
    assert!(MIN_MAPQ_UNIQ <= MAX_MAPQ);

    let alignment_probability = {
        let prob_best = best_alignment.alignment_score.exp2();
        if best_alignment_interval_size > 1 {
            // Multi-mapping
            1.0 / best_alignment_interval_size as f32
        } else if other_alignments.is_empty() {
            // Unique mapping
            1.0
        } else {
            // Pseudo-unique mapping
            let weighted_suboptimal_alignments = other_alignments
                .iter()
                // Filtering out "same" hits e.g. caused by InDels in homopolymers
                .filter(|suboptimal_alignment| {
                    suboptimal_alignment.interval != best_alignment.interval
                })
                .fold(0.0, |acc, suboptimal_alignment| {
                    suboptimal_alignment
                        .alignment_score
                        .exp2()
                        .mul_add(suboptimal_alignment.interval.size as f32, acc)
                });
            prob_best / (prob_best + weighted_suboptimal_alignments)
        }
    }
    // Guard against rounding errors
    .clamp(0.0, 1.0);

    // Produce Phred score
    let mapping_quality = (-10.0 * (1.0 - alignment_probability).log10())
        .min(f32::from(MAX_MAPQ))
        .round() as u8;

    // When we found a best-scoring hit, we search up until `AS + mm_penalty` to find suboptimal
    // hits in order to estimate the mapping quality. However, we don't ever search beyond the
    // `mismatch_bound` limit. So when we are closer than `mm_penalty` to the `mismatch_bound`,
    // the mapping quality needs to be scaled down to reflect the fact that we are less likely to
    // find suboptimal alignments within these narrowed bounds.
    if mapping_quality == MAX_MAPQ {
        let scaled_mq = f32::from(MAX_MAPQ - MIN_MAPQ_UNIQ).mul_add(
            alignment_parameters
                .mismatch_bound
                .remaining_frac_of_repr_mm(
                    best_alignment.alignment_score,
                    best_alignment.edit_operations.read_len(),
                )
                .min(1.0),
            f32::from(MIN_MAPQ_UNIQ),
        );
        return scaled_mq.round() as u8;
    }

    mapping_quality
}

/// Create and return a BAM record of either a hit or an unmapped read
#[allow(clippy::too_many_arguments)]
fn create_bam_record(
    input_record: Record,
    relative_position: Option<u64>,
    absolute_position: Option<usize>,
    hit_interval: Option<&HitInterval>,
    mapq: Option<u8>,
    tid: Option<u32>,
    strand: Option<Direction>,
    duration: Option<&Duration>,
    // Contains valid content for the `YA` tag
    alternative_hits: Option<AlternativeAlignments>,
    original_symbols: &OriginalSymbols,
) -> Result<sam::alignment::RecordBuf> {
    let mut bam_builder = sam::alignment::RecordBuf::builder();

    let (cigar, md_tag, edit_distance) = hit_interval.map_or((None, None, None), |hit_interval| {
        let (cigar, md_tag, edit_distance) = hit_interval.edit_operations.to_bam_fields(
            strand.expect("This is not expected to fail"),
            absolute_position.expect("to be set as `hit_interval` is `Some()`"),
            original_symbols,
        );
        (Some(cigar), Some(md_tag), Some(edit_distance))
    });

    // Copy flags from input record
    let mut flags = sam::alignment::record::Flags::from(input_record.bam_flags);

    // Remove flags mapAD does not know about (yet)
    flags.remove(sam::alignment::record::Flags::MATE_UNMAPPED);
    flags.remove(sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED);
    flags.remove(sam::alignment::record::Flags::PROPERLY_SEGMENTED);
    flags.remove(sam::alignment::record::Flags::SECONDARY);
    flags.remove(sam::alignment::record::Flags::SUPPLEMENTARY);

    if let Some(position) = relative_position {
        flags.remove(sam::alignment::record::Flags::UNMAPPED);
        bam_builder = bam_builder.set_alignment_start(
            usize::try_from(position + 1)
                .ok()
                .and_then(|pos| pos.try_into().ok())
                .ok_or_else(|| Error::InvalidIndex("Could not compute valid coordinate".into()))?,
        );
    } else {
        flags.insert(sam::alignment::record::Flags::UNMAPPED);
        flags.remove(sam::alignment::record::Flags::REVERSE_COMPLEMENTED);
        flags.remove(sam::alignment::record::Flags::PROPERLY_SEGMENTED);
    }

    // Flag read that maps to reverse strand
    if strand == Some(Direction::Backward) {
        flags.insert(sam::alignment::record::Flags::REVERSE_COMPLEMENTED);
    } else {
        flags.remove(sam::alignment::record::Flags::REVERSE_COMPLEMENTED);
    }

    // Set mandatory properties of the BAM record
    if let Some(read_name) = input_record.name {
        bam_builder = bam_builder.set_name(read_name.into());
    }
    bam_builder = bam_builder.set_flags(flags);

    // Add optional fields (or do not)

    if let Some(tid) = tid {
        let tid = i32::try_from(tid)
            .ok()
            .and_then(|tid| tid.try_into().ok())
            .ok_or_else(|| Error::InvalidIndex("Invalid reference map".into()))?;
        bam_builder = bam_builder.set_reference_sequence_id(tid);
    }

    // Some (not all) fields need to be reversed when the read maps to the reverse strand
    match strand {
        Some(Direction::Forward) | None => {
            bam_builder = bam_builder
                .set_sequence(sam::alignment::record_buf::Sequence::from(
                    input_record.sequence,
                ))
                .set_quality_scores(input_record.base_qualities.into());
        }
        Some(Direction::Backward) => {
            // CIGAR strings and MD tags are reversed during generation
            bam_builder = bam_builder
                .set_sequence(sam::alignment::record_buf::Sequence::from(dna::revcomp(
                    &input_record.sequence,
                )))
                .set_quality_scores(
                    input_record
                        .base_qualities
                        .iter()
                        .rev()
                        .copied()
                        .collect::<Vec<_>>()
                        .into(),
                );
        }
    }

    // CIGAR strings and MD tags are reversed during generation
    if let Some(cigar) = cigar {
        bam_builder = bam_builder.set_cigar(cigar.into());
    }

    if let Some(mapq) = mapq {
        bam_builder = bam_builder
            .set_mapping_quality(mapq.try_into().expect("This is not expected to happen"));
    }

    // Add optional tags

    // BAM auxiliary data tags that will be removed
    let tag_filter = [
        b"AS", b"MD", b"NM", b"X0", b"X1", b"XA", b"XD", b"XE", b"XF", b"XG", b"XM", b"XN", b"XO",
        b"XS", b"XT",
    ];

    // Add tags that were already present in the input
    let mut aux_data = input_record
        .bam_tags
        .into_iter()
        // Remove BWA (+ mapAD) specific auxiliary fields (avoids potential confusion)
        .filter(|(tag, _v)| !tag_filter.contains(&tag))
        .map(|(tag, value)| Ok((tag.into(), value.into())))
        .collect::<Result<Vec<_>>>()?;

    if let Some(hit_interval) = hit_interval {
        aux_data.push((
            sam::alignment::record::data::field::tag::Tag::ALIGNMENT_SCORE,
            sam::alignment::record_buf::data::field::Value::Float(hit_interval.alignment_score),
        ));
    };

    if let Some(edit_distance) = edit_distance {
        aux_data.push((
            sam::alignment::record::data::field::tag::Tag::EDIT_DISTANCE,
            sam::alignment::record_buf::data::field::Value::Int32(i32::from(edit_distance)),
        ));
    };

    // CIGAR strings and MD tags are reversed during generation
    if let Some(md_tag) = md_tag {
        aux_data.push((
            sam::alignment::record::data::field::Tag::MISMATCHED_POSITIONS,
            sam::alignment::record_buf::data::field::Value::String(BString::from(md_tag)),
        ));
    }

    // Alternative alignments (augmented BWA-style)
    if let Some(alternative_hits) = alternative_hits {
        if !alternative_hits.xa.is_empty() {
            aux_data.push((
                (*b"XA").into(),
                sam::alignment::record_buf::data::field::Value::String(alternative_hits.xa.into()),
            ));
        }
        aux_data.push((
            (*b"X0").into(),
            sam::alignment::record_buf::data::field::Value::Int32(alternative_hits.x0),
        ));
        aux_data.push((
            (*b"X1").into(),
            sam::alignment::record_buf::data::field::Value::Int32(alternative_hits.x1),
        ));
        if alternative_hits.x1 > 0 {
            aux_data.push((
                (*b"XS").into(),
                sam::alignment::record_buf::data::field::Value::Float(alternative_hits.xs),
            ));
        }
        aux_data.push((
            (*b"XT").into(),
            sam::alignment::record_buf::data::field::Value::Character(
                alternative_hits
                    .xt
                    .try_into()
                    .expect("Char is guaranteed to be ASCII"),
            ),
        ));
    }

    if let Some(duration) = duration {
        // Add the time that was needed for mapping the read
        aux_data.push((
            (*b"XD").into(),
            sam::alignment::record_buf::data::field::Value::Float(duration.as_secs_f32()),
        ));
    }

    bam_builder = bam_builder.set_data(
        aux_data
            .into_iter()
            .collect::<sam::alignment::record_buf::Data>(),
    );

    Ok(bam_builder.build())
}

/// Checks stop-criteria of stack frames before pushing them onto the stack.
/// Since push operations on heaps are costly, this should accelerate the alignment.
#[allow(clippy::too_many_arguments)]
fn check_and_push_stack_frame<MB>(
    mut stack_frame: MismatchSearchStackFrame,
    pattern: &[u8],
    edit_operation: EditOperation,
    edit_tree: &mut Tree<EditOperation>,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    intervals: &mut BinaryHeap<HitInterval>,
    mismatch_bound: &MB,
    alignment_parameters: &AlignmentParameters,
) where
    MB: MismatchBound,
{
    // This is technically redundant. However, our benchmarks suggest
    // that having this here improves mapping speed
    if let Some(best_scoring_interval) = intervals.peek() {
        if mismatch_bound.reject_iterative(
            stack_frame.alignment_score,
            best_scoring_interval.alignment_score,
        ) {
            return;
        }
    }

    if stack_frame.num_gaps_open > alignment_parameters.max_num_gaps_open {
        return;
    }

    stack_frame.edit_node_id = edit_tree
        .add_node(edit_operation, stack_frame.edit_node_id)
        .expect("We bound the length of `edit_tree` at `STACK_LIMIT` < `u32`");

    if stack_frame.backward_index < 0 && stack_frame.forward_index > (pattern.len() as i16 - 1) {
        // This route through the read graph is finished successfully, push the interval
        let edit_operations =
            extract_edit_operations(stack_frame.edit_node_id, edit_tree, pattern.len());
        intervals.push(HitInterval {
            interval: stack_frame.current_interval,
            alignment_score: stack_frame.alignment_score,
            edit_operations,
        });
        //print_debug(&stack_frame, intervals, edit_tree);
        return;
    }

    stack.push(stack_frame);
}

/// Prints backtracking trace to stderr
#[allow(dead_code)]
fn print_debug(
    stack_frame: &MismatchSearchStackFrame,
    intervals: &BinaryHeap<HitInterval>,
    edit_tree: &Tree<EditOperation>,
) {
    let best_as = intervals.peek().map_or(0.0, |v| v.alignment_score);

    eprintln!(
        "{}\t{}\t{}\t{}\t{}\t{:?}",
        stack_frame.alignment_score,
        stack_frame.alignment_score,
        stack_frame.backward_index,
        stack_frame.forward_index,
        best_as,
        edit_tree.ancestors(stack_frame.edit_node_id).next(),
    );
}

/// Finds all suffix array intervals for the current pattern
/// w.r.t. supplied alignment parameters
#[allow(clippy::too_many_arguments)]
pub fn k_mismatch_search<SDM, MB>(
    pattern: &[u8],
    base_qualities: &[u8],
    parameters: &AlignmentParameters,
    fmd_index: &RtFmdIndex,
    stack: &mut MinMaxHeap<MismatchSearchStackFrame>,
    edit_tree: &mut Tree<EditOperation>,
    sequence_difference_model: &SDM,
    mismatch_bound: &MB,
) -> BinaryHeap<HitInterval>
where
    SDM: SequenceDifferenceModel,
    MB: MismatchBound,
{
    let center_of_read = pattern.len() / 2;
    let bi_d_array = BiDArray::new(
        pattern,
        base_qualities,
        center_of_read,
        parameters,
        fmd_index,
        sequence_difference_model,
    );

    let optimal_penalties =
        compute_optimal_scores(pattern, base_qualities, sequence_difference_model);
    let mut hit_intervals: BinaryHeap<HitInterval> = BinaryHeap::new();

    let mut stack_size_limit_reported = false;
    stack.clear();
    let root_node = edit_tree.clear();

    stack.push(MismatchSearchStackFrame {
        current_interval: fmd_index.init_interval(),
        backward_index: center_of_read as i16 - 1,
        forward_index: center_of_read as i16,
        direction: Direction::Forward,
        gap_backwards: GapState::Closed,
        gap_forwards: GapState::Closed,
        num_gaps_open: 0,
        alignment_score: 0.0,
        edit_node_id: root_node,
    });

    let mut mm_scores = [0_f32; 4];

    while let Some(stack_frame) = stack.pop_max() {
        // Determine direction of progress for next iteration on this stack frame
        let (
            j,
            next_backward_index,
            next_forward_index,
            fmd_ext_interval,
            optimal_penalty,
            next_insertion_backward,
            next_insertion_forward,
            next_deletion_backward,
            next_deletion_forward,
            insertion_score,
            deletion_score,
            next_closed_gap_backward,
            next_closed_gap_forward,
            num_gaps_open,
        );
        match stack_frame.direction {
            Direction::Forward => {
                next_forward_index = stack_frame.forward_index + 1;
                next_backward_index = stack_frame.backward_index;
                j = stack_frame.forward_index;
                optimal_penalty = optimal_penalties[j as usize];
                fmd_ext_interval = stack_frame.current_interval.swapped();
                next_insertion_backward = stack_frame.gap_backwards;
                next_insertion_forward = GapState::Insertion;
                next_deletion_backward = stack_frame.gap_backwards;
                next_deletion_forward = GapState::Deletion;
                next_closed_gap_backward = stack_frame.gap_backwards;
                next_closed_gap_forward = GapState::Closed;
                insertion_score = if stack_frame.gap_forwards == GapState::Insertion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open + parameters.penalty_gap_extend
                } + stack_frame.alignment_score;
                deletion_score = if stack_frame.gap_forwards == GapState::Deletion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open + parameters.penalty_gap_extend
                } + stack_frame.alignment_score;
                for (score, &base) in mm_scores.iter_mut().zip(b"ACGT".iter().rev()) {
                    *score = sequence_difference_model.get(
                        j as usize,
                        pattern.len(),
                        dna::complement(base),
                        pattern[j as usize],
                        base_qualities[j as usize],
                    ) - optimal_penalty
                        + stack_frame.alignment_score;
                }
                num_gaps_open = if stack_frame.gap_forwards == GapState::Closed {
                    stack_frame.num_gaps_open + 1
                } else {
                    stack_frame.num_gaps_open
                };
            }
            Direction::Backward => {
                next_forward_index = stack_frame.forward_index;
                next_backward_index = stack_frame.backward_index - 1;
                j = stack_frame.backward_index;
                optimal_penalty = optimal_penalties[j as usize];
                fmd_ext_interval = stack_frame.current_interval;
                next_insertion_backward = GapState::Insertion;
                next_insertion_forward = stack_frame.gap_forwards;
                next_deletion_backward = GapState::Deletion;
                next_deletion_forward = stack_frame.gap_forwards;
                next_closed_gap_backward = GapState::Closed;
                next_closed_gap_forward = stack_frame.gap_forwards;
                insertion_score = if stack_frame.gap_backwards == GapState::Insertion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open + parameters.penalty_gap_extend
                } + stack_frame.alignment_score;
                deletion_score = if stack_frame.gap_backwards == GapState::Deletion {
                    parameters.penalty_gap_extend
                } else {
                    parameters.penalty_gap_open + parameters.penalty_gap_extend
                } + stack_frame.alignment_score;
                for (score, &base) in mm_scores.iter_mut().zip(b"ACGT".iter().rev()) {
                    *score = sequence_difference_model.get(
                        j as usize,
                        pattern.len(),
                        base,
                        pattern[j as usize],
                        base_qualities[j as usize],
                    ) - optimal_penalty
                        + stack_frame.alignment_score;
                }
                num_gaps_open = if stack_frame.gap_backwards == GapState::Closed {
                    stack_frame.num_gaps_open + 1
                } else {
                    stack_frame.num_gaps_open
                };
            }
        };

        // Calculate the lower bounds for extension
        let lower_bound = bi_d_array.get(next_backward_index, next_forward_index);

        //print_debug(&stack_frame, &hit_intervals, edit_tree);

        // Since we operate on a priority stack, we can assume that there are no
        // better scoring frames on the stack, so we are going to stop the search.
        if let Some(best_scoring_interval) = hit_intervals.peek() {
            if mismatch_bound.reject_iterative(
                stack_frame.alignment_score + lower_bound,
                best_scoring_interval.alignment_score,
            ) {
                break;
            }
        }

        //
        // Insertion in read / deletion in reference
        //
        {
            if !mismatch_bound.reject(insertion_score + lower_bound, pattern.len())
                && i16::from(parameters.gap_dist_ends) <= j.min(pattern.len() as i16 - j)
            {
                check_and_push_stack_frame(
                    MismatchSearchStackFrame {
                        backward_index: next_backward_index,
                        forward_index: next_forward_index,
                        direction: stack_frame.direction.reverse(),
                        // Mark opened gap at the corresponding end
                        gap_backwards: next_insertion_backward,
                        gap_forwards: next_insertion_forward,
                        alignment_score: insertion_score,
                        num_gaps_open,
                        ..stack_frame
                    },
                    pattern,
                    EditOperation::Insertion(j as u16),
                    edit_tree,
                    stack,
                    &mut hit_intervals,
                    mismatch_bound,
                    parameters,
                );
            }
        }

        // Bidirectional extension of the (hit) interval
        for ((c, mut interval_prime), mm_score) in fmd_index
            .extend_iter(&fmd_ext_interval)
            .zip(mm_scores.iter())
        {
            if interval_prime.size < 1 {
                continue;
            }

            // Special treatment of forward extension
            let c = match stack_frame.direction {
                Direction::Forward => {
                    interval_prime = interval_prime.swapped();
                    dna::complement(fmd_index.get_rev(c))
                }
                Direction::Backward => fmd_index.get_rev(c),
            };

            //
            // Deletion in read / insertion in reference
            //
            {
                if !mismatch_bound.reject(deletion_score + lower_bound, pattern.len())
                    && i16::from(parameters.gap_dist_ends) <= j.min(pattern.len() as i16 - j)
                {
                    check_and_push_stack_frame(
                        MismatchSearchStackFrame {
                            current_interval: interval_prime,
                            // Mark open gap at the corresponding end
                            gap_backwards: next_deletion_backward,
                            gap_forwards: next_deletion_forward,
                            alignment_score: deletion_score,
                            num_gaps_open,
                            ..stack_frame
                        },
                        pattern,
                        EditOperation::Deletion(j as u16, c),
                        edit_tree,
                        stack,
                        &mut hit_intervals,
                        mismatch_bound,
                        parameters,
                    );
                }
            }

            //
            // Match/mismatch
            //
            {
                if !mismatch_bound.reject(mm_score + lower_bound, pattern.len()) {
                    check_and_push_stack_frame(
                        MismatchSearchStackFrame {
                            current_interval: interval_prime,
                            backward_index: next_backward_index,
                            forward_index: next_forward_index,
                            direction: stack_frame.direction.reverse(),
                            // Mark closed gap at the corresponding end
                            gap_backwards: next_closed_gap_backward,
                            gap_forwards: next_closed_gap_forward,
                            alignment_score: *mm_score,
                            ..stack_frame
                        },
                        pattern,
                        if c == pattern[j as usize] {
                            EditOperation::Match(j as u16)
                        } else {
                            EditOperation::Mismatch(j as u16, c)
                        },
                        edit_tree,
                        stack,
                        &mut hit_intervals,
                        mismatch_bound,
                        parameters,
                    );
                }
            }
        }

        // Only search until we've found a multi-hit

        // The best scoring hit is always guaranteed to be at the top of the priority queue,
        // but the second best hit is not guaranteed to be second of the priority queue at all
        // times. This is because with one source stack frame, we can find and store up to 9 hits!
        // Hence, the second best hit is only guaranteed to be found with the next best
        // second-to-last source stack frame.
        if (hit_intervals.len() > 9)
            || hit_intervals
                .peek()
                .filter(|best_hit| best_hit.interval.size > 1)
                .is_some()
        {
            return hit_intervals;
        }

        // Limit stack size
        if stack.len() > STACK_LIMIT as usize || edit_tree.len() > EDIT_TREE_LIMIT as usize {
            if parameters.stack_limit_abort {
                return hit_intervals;
            }
            if !stack_size_limit_reported {
                trace!(
                    "Stack size limit exceeded (read length: {} bp). Remove highly penalized partial alignments from stack (stack size: {}, edit tree size: {}).",
                    pattern.len(),
                    stack.len(),
                    edit_tree.len(),
                );
                stack_size_limit_reported = true;
            }

            for _ in 0..(stack.len() as isize - STACK_LIMIT as isize)
                .max(edit_tree.len() as isize - EDIT_TREE_LIMIT as isize)
            {
                // The stack should never be empty at this point, so we boldly ignore the
                // `None` case here
                if let Some(min_stack_frame) = stack.pop_min() {
                    edit_tree.remove(min_stack_frame.edit_node_id);
                }
            }
        }
    }
    hit_intervals
}

#[cfg(test)]
pub mod tests {
    use std::collections::BTreeMap;

    use assert_approx_eq::assert_approx_eq;
    use bio::alphabets;
    use noodles::sam;

    use super::*;
    use crate::{
        index::DNA_UPPERCASE_ALPHABET,
        map::{mismatch_bounds::*, sequence_difference_models::*, *},
        utils::*,
    };

    #[test]
    fn test_inexact_search() {
        let difference_model = TestDifferenceModel {
            deam_score: -0.5,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "ACGTACGTACGTACGT".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GTTC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score, vec![-1.0]);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![2, 6, 10, 19, 23, 27]);
    }

    #[test]
    fn test_reverse_strand_search() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -20.0,
            penalty_gap_extend: -10.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "GAAAAG".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // 1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![8]);
    }

    #[test]
    fn test_gapped_alignment() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -3.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "TAT".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -2.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0, 2, 5]);
    }

    #[test]
    fn test_gapped_alignment_read_end() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -6.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "AAAAAGGGGAAAAA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        // Gap in the middle of the read (allowed)
        let pattern = "AAAAAAAAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        assert_eq!(positions, vec![0]);

        // Gap near read end (not allowed)
        let pattern = "AGGGAAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        assert_eq!(positions, vec![]);
    }

    #[test]
    fn test_gap_open_limit() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -6.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
            max_num_gaps_open: 1,
        };

        let ref_seq = "CTAGCCAGCGATTTACATGCTCTCGGAATATCGACATGTA"
            .as_bytes()
            .to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        // One gap (allowed)
        let pattern = "CTAGCCAGCGAACATGCTCTCGGAATATCGACATGTA"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        assert!(positions.contains(&0));

        // Two gap (not allowed)
        let pattern = "CTAGCCAGCGATTACATGCTCTCGGAATTCGACATGTA"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&suffix_array))
            .collect();
        assert!(positions.is_empty());
    }

    #[test]
    fn test_vindija_pwm_alignment() {
        let difference_model = VindijaPwm::new();
        let mmb = TestBound {
            threshold: -30.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            // Disable gaps
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -200.0,
            penalty_gap_extend: -100.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "CCCCCC".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTCCCT".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score[0], -4.641691);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&sar))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0]);

        let pattern = "CCCCCC".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_eq!(alignment_score[0], 0.0);

        let mut positions: Vec<usize> = intervals
            .into_iter()
            .flat_map(|f| f.interval.forward().occ(&sar))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0]);

        //
        // Test "normal" mismatch
        //

        let ref_seq = "AAAAAA".as_bytes().to_owned(); // revcomp = "ATA"

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "AAGAAA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -30.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_score: Vec<f32> = intervals.iter().map(|f| f.alignment_score).collect();
        assert_approx_eq!(alignment_score[0], -10.965062);
    }

    #[test]
    fn test_ord_impl() {
        let mut edit_tree: Tree<usize> = Tree::new();
        let root_id = edit_tree.clear();

        let map_params_large = MismatchSearchStackFrame {
            alignment_score: -5.0,
            current_interval: RtBiInterval {
                lower: 5,
                lower_rev: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            direction: Direction::Backward,
            gap_forwards: GapState::Closed,
            gap_backwards: GapState::Closed,
            num_gaps_open: 0,
            edit_node_id: root_id,
        };
        let map_params_small = MismatchSearchStackFrame {
            alignment_score: -20.0,
            current_interval: RtBiInterval {
                lower: 5,
                lower_rev: 5,
                size: 5,
            },
            backward_index: 5,
            forward_index: 5,
            direction: Direction::Backward,
            gap_forwards: GapState::Closed,
            gap_backwards: GapState::Closed,
            num_gaps_open: 0,
            edit_node_id: root_id,
        };

        assert!(map_params_large > map_params_small);
        assert!(map_params_large >= map_params_small);
        assert_ne!(map_params_large, map_params_small);
    }

    #[test]
    fn test_corner_cases() {
        let difference_model = VindijaPwm::new();
        let repr_mm_penalty = difference_model.get_representative_mismatch_penalty();
        let mmb = Discrete::new(0.01, 0.02, repr_mm_penalty);

        let parameters = AlignmentParameters {
            penalty_gap_open: 3.0 * difference_model.get_representative_mismatch_penalty(),
            penalty_gap_extend: 0.6 * difference_model.get_representative_mismatch_penalty(),
            difference_model: difference_model.clone().into(),
            chunk_size: 1,
            mismatch_bound: mmb.clone().into(),
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        // "correct" "AAAAAAAAAAAAAAAAAAAA" (20x 'A') "incorrect"
        let ref_seq = "GTTGTATTTTTAGTAGAGACAGGGTTTCATCATGTTGGCCAGAAAAAAAAAAAAAAAAAAAATTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GTTGTATTTTTAGTAGAGACAGGCTTTCATCATGTTGGCCAG"
            .as_bytes()
            .to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let alignment_scores = intervals
            .iter()
            .map(|f| f.alignment_score)
            .collect::<Vec<_>>();
        assert_eq!(alignment_scores, vec![-10.936638, -39.474224, -10.965062]);

        let mut positions: Vec<usize> = intervals
            .iter()
            .flat_map(|f| f.interval.forward().occ(&sar))
            .collect();
        positions.sort();
        assert_eq!(positions, vec![0, 62, 63]);

        assert_eq!(
            intervals.peek().unwrap().interval.forward().occ(&sar),
            vec![0]
        );
    }

    #[test]
    fn test_cigar_indels() {
        let difference_model = TestDifferenceModel {
            deam_score: -10.0,
            mm_score: -10.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -4.0,
            representative_mm_bound: -10.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        //
        // Deletion
        //
        let ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(
            cigar,
            vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    4
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Deletion,
                    1
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    2
                ),
            ]
        );

        //
        // 2-base deletion
        //
        let ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(best_hit.alignment_score, -4.0);
        assert_eq!(
            cigar,
            vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    3
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Deletion,
                    2
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    3
                ),
            ]
        );

        //
        // Insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(best_hit.alignment_score, -3.0);
        assert_eq!(
            cigar,
            vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    5
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Insertion,
                    1
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    2
                ),
            ]
        );

        //
        // 2-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(best_hit.alignment_score, -4.0);
        assert_eq!(
            cigar,
            vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    5
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Insertion,
                    2
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    2
                ),
            ]
        );

        //
        // 3-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGTGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mmb = TestBound {
            threshold: -5.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -4.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (cigar, _, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(best_hit.alignment_score, -5.0);
        assert_eq!(
            cigar,
            vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    5
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Insertion,
                    3
                ),
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    2
                ),
            ]
        );
    }

    #[test]
    fn test_md_tag() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -2.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -2.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        //
        // Mutation
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTATA".as_bytes().to_owned();
        let base_qualities = vec![40; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(md_tag, "5C1".as_bytes());

        //
        // Deletion
        //
        let ref_seq = "GATTAGCA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "ATTACA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mmb = TestBound {
            threshold: -4.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -2.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(md_tag, "4^G2".as_bytes());

        //
        // 2-base deletion
        //
        let ref_seq = "GATTACAG".as_bytes().to_owned(); // CTGTAATC

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATCAG".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(md_tag, "3^TA3".as_bytes());

        //
        // Insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(md_tag, "7".as_bytes());

        //
        // 2-base insertion
        //
        let ref_seq = "GATTACA".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _sar) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "GATTAGGCA".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -3.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );
        let best_hit = intervals.pop().unwrap();
        let (_, md_tag, _) = best_hit.edit_operations.to_bam_fields(
            Direction::Forward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );

        assert_eq!(md_tag, "7".as_bytes());
    }

    #[test]
    fn test_reverse_strand_search_2() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: 0.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "AAAGCGTTTGCG".as_bytes().to_owned();

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TTT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let mut intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // 0.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let positions = if let Some(best_alignment) = intervals.pop() {
            best_alignment
                .interval
                .forward()
                .occ(&suffix_array)
                .iter()
                .filter(|&&position| position < suffix_array.len() / 2)
                .map(|&position| (position, Direction::Forward))
                .chain(
                    best_alignment
                        .interval
                        .revcomp()
                        .occ(&suffix_array)
                        .iter()
                        .filter(|&&position| position < suffix_array.len() / 2)
                        .map(|&position| (position, Direction::Backward)),
                )
                .collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        assert_eq!(
            positions,
            vec![(6, Direction::Forward), (0, Direction::Backward),]
        );
    }

    #[test]
    fn test_edit_operations_reverse_strand() {
        let difference_model = TestDifferenceModel {
            deam_score: -1.0,
            mm_score: -1.0,
            match_score: 0.0,
        };
        let mmb = TestBound {
            threshold: -1.0,
            representative_mm_bound: -1.0,
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mmb.clone().into(),
            penalty_gap_open: -3.0,
            penalty_gap_extend: -1.0,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let ref_seq = "GATTACA".as_bytes().to_owned(); // revcomp = "TGTAATC"

        // Reference
        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, suffix_array) = build_auxiliary_structures(ref_seq, alphabet);

        let pattern = "TAGT".as_bytes().to_owned();
        let base_qualities = vec![0; pattern.len()];

        let mut stack_buf = MinMaxHeap::new();
        let mut tree_buf = Tree::new();
        let intervals = k_mismatch_search(
            &pattern,
            &base_qualities,
            // -1.0,
            &parameters,
            &fmd_index,
            &mut stack_buf,
            &mut tree_buf,
            &difference_model,
            &mmb,
        );

        let best_alignment = intervals.peek().unwrap();

        let positions = best_alignment
            .interval
            .forward()
            .occ(&suffix_array)
            .iter()
            .filter(|&&position| position < suffix_array.len() / 2)
            .map(|&position| (position, Direction::Forward))
            .chain(
                best_alignment
                    .interval
                    .revcomp()
                    .occ(&suffix_array)
                    .iter()
                    .filter(|&&position| position < suffix_array.len() / 2)
                    .map(|&position| (position, Direction::Backward)),
            )
            .collect::<Vec<_>>();
        assert_eq!(positions, vec![(1, Direction::Backward),]);

        let (_cigar, md, edop) = best_alignment.edit_operations.to_bam_fields(
            Direction::Backward,
            0,
            &OriginalSymbols::new(BTreeMap::new()),
        );
        assert_eq!(md, b"1T2");

        assert_eq!(edop, 1);
    }

    #[test]
    fn test_n() {
        let ref_seq = "GATTACAGATTACAGATTACA".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound = TestBound {
            threshold: -14.0,
            representative_mm_bound: difference_model.get_representative_mismatch_penalty(),
        };

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mismatch_bound.clone().into(),
            penalty_gap_open: 0.001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 0,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, _) = build_auxiliary_structures(ref_seq, alphabet);

        {
            let pattern = "NNNNNNNNNN".as_bytes().to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        {
            let pattern = "AGATNACAG".as_bytes().to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }
    }

    #[test]
    fn test_bench() {
        let ref_seq = "AACATCTTTTTCGGCCGTTTCATAAATACTCTTTATTCACAATTAAGGTCTCTCTGTGCTTCTACGATGGTTGGCTGGCCTGCTACTATGTAGCCCAGTT\
TTGTACTTGAAACAGATGAGTCTCTTCACACTCATACGCGTTTCAGCAAGTACCAGTGATAGCGGTGGTGCTCCTCAATATCCGCCGCTCGCAAGCGACC\
GAGAAAAGCACGCATCGAGAACCTTGTCTTACACAACTTTCAATGCTAGCGCCTAAGTTCAACCGCTAATATAACAGCTGAGTCCCCATCGCCCGAGGGT\
CGCATAGTTACTCCGTGAGGTGTATCAGAAATCCGCGATTTCGTAACTTCTACATCCCGTCCTCTCAGACCACCATTTTCCGTCATTACCAGCAAATATG\
ATCCGCCGATCGGCCGGACATCTGCCAAACCTTGCTGTCCTCCAATAAGTTTAGGTTAAGGGCCCAGTCAAGAAAAGACCCCAGTTCAAATGACGTACAA\
GGCAGCACCCCGTTTATACACCGGGGCATGCTTAGCTGCTATCGGACTTGCGCGGATAAAAATAGCCTAACTGCTAAAAGCATGTCCGCCCACGATGGAC\
TGGGGATTAACAGATAAGTTTAACATCTTAAACGTGAGTTGCGACAACAGCCAGGCTGACCGTACTGCCCGGACGTGCGGGGTCTGCTCATTGTGAATTC\
AGAGATGCCTCTATATTAGATTGGAGGGTAACTCAATAATCATGAGCGCAAGCTCACGAGAAGTACTGCTACTCCGCATATGTTGCTGAATATCTTCTGC\
CTGCGTGGCGGCCAGCATGTGATTACAATGACTCACTAGGGGGGGTAGTAGAGTCCGATTCGGTTGACTACTCGGGGTAACTAGGCATACAGTGGGCTCT\
AAGTCTGTGCGTCAAGCCCAGGGCAGGTCGTCGCTCATTAAGTTAGCCCACTTCATTAAGTATCACAGGGGCATCCTAGGTCCGTCTTAGTGTCCAACCA\
GCACCTCAATGTAATGTGGTGAGGGTAGTCACCAGTAGCGGAGTCCATCCTCTGGCCGCGTCGAGTCTTTCACGTCTATTCATGACGTATAGTCCCATTG\
GTGCGTCATGTACAGGCTTCTAGGACCACATTCTCCTTGGGCGGTTATCAAGCAAAATGCAAAGGTACCTCGGCATTCCGCCGATATTGTCTCTGAGATG\
TACGACGAGCCACACACTACGCTCGTTCTTCTGAGCCAATGGTTCGTATGGGATGCACTATAGTCGGACTAGACCCTTGCCATACTACATGCTGAAGTTT\
CAGGGATAACAATCGGTGGGCATAGACAACGCTAACGTCCGCGTTAAAGATTTAGTAACGTACGATACGCTGATTGCTGCAGGACAGAATAGACTTCAAT\
CGCGCTACCCTTCTAGCGCTGGGTAGTTTCGCGGACCGTATGCCAGGTAACAAGCACAGTAAACGCACTATTGGGTAATACGGATCTCAAAGGCTCATGC\
TCGGTGAAGGTGCCCGTAAGTGCACTTTCCTAACATAGCCGGCGTTATTAAGAGTTCCATCGGACAGCGCCCCACACCAAATCACTCAATTTCTCCTGTT\
CACCTTGCGCGTCATAATCTAACTGTGGCGAAGAACTAGGAACCTGAGAATAACCATTTCGATCGGTGATCCGTCAGATCGCTATTAGGAGCCAGAGACG\
TCTTCTACGGGGTCACCGTTCGTACCATGCTACCAGTCTTCGTCCCCTATGAGCATTTAACACAATCGTACTGTAGGTTGTCACACAACCGCTTAACGTT\
GTATCAGCTTTGCATACGCGAAGCAACGATGAGAACCGCTAAGGCTTCAATCTCCCGTTAAAGACCCTGAACTGGTCCGTGGCGGCGTGTTTTATCAGTG\
TACGGGATTCTTCGTGTACTTACTTGATAACGGCTCATACCAAGGCACGTGATTCTCTCGCAGAGGTCTCCCCCTTGGAAGTTATAGGCGCGCCGGTTTC\
CATTAGCATTGGACCACGACCGGGGGTGCTCCTATGTCCGACTCTCCATAGCGAGGAACATTGCTTCAAGCCCTTGGACGCTAGTTATTCATGCCTTGCG\
ATCTCCTAGCTTACCAAGGTTGTCGAACATCCGTTGAGAAAGTCATAAGACGTCAATAAACGTAACGGGGTTCTTTCTGCACTTGCGCTCAGACCTGCAT\
TTTAACCGTACAGGGAGCATCATACATCCCGTAACCACAGATATATCACTAACCAGGAATGACAAAGTCGGAAAGTAAGTCCCTACCGTCCGGAACCAGC\
GACGACCCAGGAAACTCCAGACTAGAAAGAGTACAGCAACCATGATATAAATAGTGCGAGCTGACAATGGTTTCGGCAGCACTGAATTTTTCGGGCTCAC\
AACAGGCGACAACCGGCCTTCGGAGAGTTTAATGCGCCATTCAGTTCGTCATCACGCTGGTACCGTCACCCTATGGGAGCTGCCAACTGTCACGATAGAA\
GATTATTTGTCACGGGACCGATTGCTTAAGACAATAACCAGGCTGGGACTACGTCATACGTCTCCAACTAACAAAAGTACATCCCCTATGAGATCAAAGG\
TACCGTATAAGAGAGTATGTCTAGGGCCTTGCATGGTTTTTACTTAACGCCGCTCCCGGAGCTTGCCCGCCATCCGGTCGAGAAAGCGTACCTCAAGTTT\
TGACCCTAACCAACACCCTAGAAGTTATTTGCTTCTAGGTTTCATACTTAAACTGCAATTGACACCGTACAAGCAGTCCGACCACTTTGGAAGCGTCGTG\
GATGCCCGGTTATACGGGCGATCGAGCTCTAGTCCTGACTCACGGCCGTGTATAACTTGATTCACGCCCCTTGCGCCGAACGCCCTATGGATCAAGGCGC\
CATCGAAGAGGTATTGTTACCATCAGAGGTTCGTGTGTCCTTTTAGCGCAACATTGATTACCGGTGCCTCCGTAAGTGGGGCCGACCAAAGTGCCGCCAG\
ATACTACCCCGCACGAAAATAGAAATGCGTTTGCGACCCTAACACTGAGCGTGCTGCCACTTGAAATTAATACTCTGAGGGCGCGGATCACTACCGGCGA\
GGCAAGATTACAATACGTACACTCCGTTAATACCATAACTCTTGGGTGTAGGAACTCCGGTAATACACCCGCAGGACGTTATTGGCCGTCGCCAGCAGCT\
TAGTAACAGATATAAAATCTAGGGACTTGTCGGGCAATGTCGGACGAATTCCGTTGTACCGTCAACTCCCTATACGCAAACTTGGTAACGTTTACTAGGT\
AGACGCGTTCAATGAACGGGTCTACTACTTCACTGGCGCCTCGTTGTGGCGTACTGCACCCGTAATACCCTATACCCGGCGGTATTCATCAGTCGCTCTG\
TCCACAGCAGCCAGGGAGGCTGTTTTCCTCGGGGCTACGGGGTCAACCCGATCCCATGGGTTTGCGTAGTACAATGAAAGGTCTCCTATCTTTCTTTTTG\
CGTACGGCTGCAGGGAAGACGCGGACATCATAGGCAACAGGGCCAGAGATATCCCCGGCTTCGAGCACTAGGGTGCTATTCCCCCAACTACGGCGCCAGG\
ACACGTGCCACTTGTTGCCACCTAACATCACTGTTGACCGTCAGATTGAAGAGTACTCAATGTCCACGATTTGCGCTGTGCCAGCGGCATCAGAACGGGA\
ATACGCAGGCATGCGTGAGAGCCAAATACCGCCCTCTGTCGGCAATCAGGTCGTAAAAGTTGAGATGTGAGGGTCTTAAGTTTGCGCCTTATGGTCTGGA\
CTCCGGTAGAACTCGGACGTTTAAGGTAACTTAGGAGGCCCTCACCTTCGAGACTGGCCGAGCAGTCAAAATCGGTCGGGCCGACGTCTTATAACTTGGC\
AAAGGTGAGGGGCACGTCTGTAATATGCAGCTGACATGTGGCTAAGAGTTGTTAAAATGTGTCGTTAAGTGTTAATACCTCACATTAGGTTTTAACGACC\
ACCTTGTTAGGTTTAAATCGGCAACTTGGCATGGACTATGTACTGGAAATCGACATAACGTGGTGAAGACTCCATACCTCGACGCATACGACCCCTCAAT\
ATCCCGCATCGTGTCCATTTATCGCTAAATAGGAAAGTGACCGACGCATACGACCGCCAGGCGCACACAAAGCTGAATCCGGCCAATAAAATAGTCTCAT\
GGTCCCGCCCGAGCCCAGTCTCGGCACCCCTAGACCCGAGTGCCGATCTCGTAGGGACACACCCTGATCGATTTTTTGCTGCATGATCTTTAAAAGGCAG\
ATGGGGGCTTAACAGTCACTACACTGACACATGAGACCTCGGAAGCGGAATCACTGCGTTCCAAACTCTCAGAGACCACGGCAAAGGAGATTTCCGGGGC\
GTAGGCCTTTTCGGGATACACTACGTATGGTTCAGCGCGACAGCCCGCATATACGGTGCATTTATATCAATCTTTAGGGAGTAACCCATCAGGTGAACAT\
GTCCAGCCCCCGGGGTCGCGAGCAATTAGCCGTGATCGGATAGGTCCAAACATGACTGATGTAGGGCCACCGGTGGAACCAGCCCGAGCAATGATTCGCC\
CCGGTGACTCCTACTCAAGCGGTCGGTGAGAGTATTGACCTTGTGACCAACGACGGGAGCAGCAAGCTGCGTGTTCCTTACCAACTCCGGGATATTCCTT\
TGATGTTCTCTAAGACAGACGTTGAGCATCAGCCCATATTAAGTCGAACGTGGTCCACTATCGATAAGTCTGCACAATTGGACGCTAGCAGGAATGGTGT\
GAGGACCCGAATCAGGCCCTCGAGTGGGACCTGGCGATCACTGCCGTCTCCGTTGCCGGTCGTTGTCCAACCGGTAAACCGGACAGCCCACATGTGCACG\
GTAAATGTTCCGTGCATAACATAATGCGAATGACCGCCTTCTTATGTCCGTTCCGAGGCCGGTCTTCTCACCCTTTTAGCATGTGGCAAGCCGGCACCGG\
AAGTGGACGACCGTCGCGTCTGTTCGTTGAGCGTGTATCTGGCGACTTAAACGATCCTACTAATCGCGTACGCTTGGTGCCTATAGAACCCGGGGAACTT\
GGCACGTCATGCCGGATACTTAATAGCCTGGCTCGTGTCGTGAGAACAATCTCCTGTGAAGTTACGTCTGAGTTTTAAAAGGCATCGTCGATCCCCTTCC\
CCTGCGCGCGGAGACGACTTAGCGACGTCGTGCACGACGAAAGCAGCCCTATTATAAGTCGTCCAGAAAGCCAGTATCTGATGTAAGCCAGACGCTAAGT\
TCAACCCAGTACATTTAAGACTAACCCATGTAAGCAAGTTATGTTGATCGCCCTCACAGCTATTCCCCGGCGGTGGACTCCATTCCATGCAGCTGATATA\
CAGATACATCACGCGTTAATTATGGTCTAGGTGTGTAAGCCGTTTCTAGGTTAGGCAGCAGTTAACGGCACTAGAGTATTGAAACACGGGAAACTAGGGG\
GCATCAATTTTCAGATGTCGTACCACAATTGCGAAATGGAACTGCAATTGTATCAAATGGCGCAATAGATATGTAGTGTCGTCCTGAGAATTCATTGGCA\
GTAAATCAACCTACTGAGCGGAGGATGTATACGCCAAATACATACATTAAAGTATCAATTATGCCCGGCCGAGTCTGGACAACCCGAGAAGCGTGATTTT\
GGGGAGTGGGGCGACGAGTCAGTAGCGTTTGGGCAGAACATGGCCTCCCTAGAGCCGCGACCGCGTGCTGGTACGGCACTCAAACTACATTTCTTTTGAG\
CGCGTAATATGAGAATTACACGGGTGCTAGATATGTAAGATACCACCAACTTCCCGCGTCCAAGGAACGGGAGACCACCTCCCGGGCCGTGCTGAGCATT\
GCCCGTTCAAGGCACCAGCAATTCTGGTGCCATAGTCGTAGAGCAGCCAACAACGGTGAAGAAGCTTACCGAAATTACAATGCGAGTGTGGCCTGCCGGC\
ACCTAAGTGTCAGCGGCCGGGGGGAGAGGATATTTTGCCGCATCCCGTCCCAAACGCGTGTTACGATCTACTGCAGATATTTGATTAGCGGTATCTTTTG\
CCCCTCCCCGCTTCAAGTTTTCATCCTAAGGGTGGGACACTTACAAGCTACAGGCCCCTTGGTAACAGGCTGTACACATCACACGTTGGAAAGCGTAATT\
ATAAAGGGTGTAGCCCTTATCTACTACTTCCAGACCGGAGTTTAACAGAGACCTGAAAGGACACTGCATTTACGAAGGGATGCCCGGGCCAACGAGAATG\
CGAACGGGTTATCGACTTAACGCGTAGTCTTCCGCTGCGGTAATTGTGTTGGATCTCCAGCTATACAGGCTGGAGCTCACTTGCCAGTCGGTAGGGTGAG\
CGGCCTTTATAGAGCGTGAGTTAGCATTGCCATTGAAATGAGATTGCTTGCCCAATACAAGTCTGTCTTCCAAGTAATAGGGCATTTGAGGCTTGAGGAA\
GCCACGCTCAGCTCTGGATGAAAAAAGTAACACGCTTCTATTAAGAAAAAGGTGAAGAAACTCCGCCTCGCAGTAGGTATCCAGCCGAACAGGGCAATCT\
CTGAGGGGACTTGTTAGACCGTCTCTGGGACTCCCATTGGCAGCAACAAAGGGATGGATCACATCTGAAGCGTCGCGGCTATTAAGAGTCTCAGCCCGGG\
CTATGCCGTATGTTTACCAAACCTCTGGTACCCATATTTCGTAGTGATGTCGCATAAAGCGCTCCTAGCAAGTTCGATCACACGTGACGAGGAGGCACGT\
TCCGCGTCATACAGACCGGGGGTGTTACTTCATTGGTACGACGGGGCACGACTCACCGAATGGGCCTGGTATTCTTCGGACGGTAGGTCCATCTACCATG\
CTCGCCACGTTTACACGGAACTTCGGGCAGCGACTCGCTAAGGGCTGCGTTAGAAACTACCTTGCAATCGGCCCAATGACTGGCAGAAATGAATCAGAAA\
CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC\
TAACCGATTCTTAACTGCCAGTCCTCCGAAATTACCACCTGTTGGGGGGCATGATTCAGAATATGAAGTAGGATGGCCATGCGAAGAGGTGGCGCTTCAT\
GCTATCATTAGTACGGGATCTATAATCGAACACGATTTGTAATAGGGCGGGAGTCAAGCGCATATCCCACGTGGTTCCAACCGCTAGCCGCCTAACAAGT\
CTGTGATGACCAATATTGGAACCTTAGTACGTGGGGCCTTGAAAATCATTTGTTCGCTGAATAATTTACTGAAGATACCAGGTTTGTCATGATCAAGGTA\
GACCCTTGCCGACTCAATCCAATCTCCGGACCTAAAGGTGGCGATCTTAAGTCACCCACTTCGGCCTTGTTAATCGCGGCGCCCTGCATCCTTGGCTTTC\
AAGGAGGGGACATAAATGATTCGCTCGCACAAACGCAAGCTGGGCGTCTCAGGATTGCCGTTATTTGCAGCTTTAGCGGTACCATGCGGCGCCTAAATTT\
AACGAACTCTCGCCGGGGTCCCTCGTTCAATGTGACTATAGCTGGCTAAGGACGTTATCAGGCATCTCCAGGTCGTATCTCCTGGCTTGCTCTATTATGA\
CCTCCCGAGTTCATTACTGCCACCAAGCGTGTTATGCCAAGTAAGGACCAACGGGCCGCTTTAATGTACTCATATCGACGCGGTCTCTATGCTATCCACA\
ACGGCATTATTAGTAATCTACAACTATGCAATGAGTTGACCTTTGGTCCACGGACCAATAACCAGGGGGCTCTATACTGGCTCCACCCATCTTTTGCATG\
TTGGCCCAAATAATGCGCCCTTCCTAACTGGATACTTCTCCTTCCTTTGCCATACTGTAAACACGACTCGGCGGAGACGAGACATACCCAAGGTGGCCCA\
ATCCACAATTAACTAAAACCCAACCGCAGGTTAGCCGCCCTATATTGTTTCGTATCGCTTTAGTGTTGGAGGGTTATATCGCTGAACACATTATTATTGC\
GAACGTTCGGTTTCCTTCCGGTCGAAAAATAGTTACTTCATATGGGCCCCACGACCGATGGAGACCCACTGATACATTCTCGACTTGATCGGTCCGGGGT\
AAGGACAACAATGTGGCCAGGCCTTAATGCGATGGGACGATCATAGGGGTAGTATTGAGCAAGTTATACTTACTTACTTGCGGTCATTCGGGTTAAAGAG\
CTTTTGCTTTTGTGAGGCAACTGTTACGTGTATCGGACATGCATCCGATGGTTGCACCAATTTTAAAATTTGGACCTAATTGTGGGAGGCGGCTGGTTAT\
TCATGGAGCCCCTCTATATGGCACGTCGAACTCGTTAGTGCTTCGTCCGGTCCTAGACCAGGGAATACCGGTGAGCAAAGATACCCCACCCGATACGCCT\
CACACATCTCTTTCCGACTGTAACGAGCTGAAGCTAGAGTTCAGGTAAAGTCGGGCTAAGGACGTACATTAAAGCGGGGTACAATTCTGCACGCACCACC\
GTCTCAGATAGGCCGCCCAATCCAAGTTCGGGAAGTGATGGGATGACTACGGCTGCAATACTGTCTGGGATCTCTAACGCGATCACCCTCTCCATCCCTT\
AGGAAAGTCATCTTCTTTCCTAAACTGCCGCTGCTTATTGGACAACCGGTGACCGACAGCGAAACCATTGTTTGTGTAGTCTAGTAGACGACCAGCTGTC\
TGGGGTAATAACCAAAATCCGAATTCCACTATCCGGGGAGTCATAGCATACTGGGGTTAATCGAGCGCAGTCGAACGCAAAACTACACACCGAGTGAGGA\
CTCTCTCCAACGACCCTAATAGACTTCTGGAGCATCTTGCGAATAACGGACGGTATGTCGATGCTACCCCGCCTTCTGCTAGCCGGCATCGGCAACACGA\
TGAGCGTTCACGAAATTATTGAGAGAACAATTGCCTTAATACTAATTCGGGATATAGTTCAAAATGGACTGATATTCTCCACTATCTCTTACTAGACCCC\
ATTTCGGATAGAAGGCCATCGGAATAAAGCTGCACGAGGATCCGTCCCGACACAATCGGAGTAGAGGCTGGACATTCGATACAGGCCGTGGGCTAGCGAC\
CTCGCTCGAGCTCAAACGTTTCAGGCTGTCGGAACCCAGACGCCCCCCGAGGCCTCTGTGGTGGCAATATCCACAATAGTGCTTAGCGTCGATCCTCTGG\
CCATTCCCGTAACCCATCATAGTAACTGTGCGGATTCTCATCCTATGGGCAATAAAGCTCCAGGTCGTTCCTTATTTCCTCTCTAACTCCCCGTTGTGTC\
GTCGTTGGACGGCCTTCGTGCCGTGTTGGAATGTAAGATAAATCACGGGGCCATAATGCATTAATATTGGACGAGCTGTTCGCAAACACGTGAAGTAAAT\
AGAATTGTAATCGAAATTGTCAAGGACTATCTGATTACCAAACTACCGACTCCCCTTCACTGCGGGCAGCTCAGGTTGTCTGTACCGGGATGATGCAAGG\
CGTTCCCCATTCAGCTTGGAGGTGGACGCAGGTTGCCAATCGCCACTCGCGCTCTAGACGTTGCATCAGGTAGGGGGGGGCTTGTTTCATAACCTAGGGG\
TTCGTGAGGGTGGATTGGAATACTGACTATCACCTAGCAATGTTTTGCAGTATTATTTTAAGATCGTGTGCGGGGTGGGGAGCTATGTTTAAATAGCCTC\
GGCTCAACTTGCCAGCGTCCCATGAACTGTACAACCTAATCGTCCCTATCGTTACTGACCGACTGCGATAGAACGCCACCTTCATAGCGTCGCCGGACAA\
GCCTGTATGCAACCCATGAGTTTCCTTCGACTAGATCCAAACTCGAGGAGGTCATGGCGAGTCAAATTGTATATCTAGCGCCCACCTGATACCTAGGTTC\
".as_bytes().to_owned();

        let difference_model = SimpleAncientDnaModel::new(
            LibraryPrep::SingleStranded {
                five_prime_overhang: 0.475,
                three_prime_overhang: 0.475,
            },
            0.001,
            0.9,
            0.02 / 3.0,
            false,
        );

        let representative_mismatch_penalty =
            difference_model.get_representative_mismatch_penalty();

        let mismatch_bound = Discrete::new(0.04, 0.02, representative_mismatch_penalty);

        let parameters = AlignmentParameters {
            difference_model: difference_model.clone().into(),
            mismatch_bound: mismatch_bound.clone().into(),
            penalty_gap_open: 0.00001_f32.log2(),
            penalty_gap_extend: representative_mismatch_penalty,
            chunk_size: 1,
            gap_dist_ends: 5,
            stack_limit_abort: false,
            max_num_gaps_open: 2,
        };

        let alphabet = alphabets::Alphabet::new(DNA_UPPERCASE_ALPHABET);
        let (fmd_index, __sar) = build_auxiliary_structures(ref_seq, alphabet);

        // bench_exogenous_read
        {
            let pattern = "GATATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTAC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        // bench_exogenous_read_deam
        {
            // 3' terminal C->T
            let pattern = "TTTATCTCGGCTGACAAACCAACAAAAAGTATCGGAACATCGCGGCGGCGTAGATGAATCTTAACCACACTCGACAGCTGTGCTTCTATACTAGCATTTT"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 0);
        }

        // bench_endogenous_read_perfect
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGACCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_1_mm_center
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTAGGGTGAGCGTAAGATTTAGGTCCATCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_2_mm_center
        {
            let pattern = "CCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_1_deam
        {
            let pattern = "TCCGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGCC"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }

        // bench_endogenous_read_4_deam
        {
            let pattern = "TCTGACAGCTATAAGGAGCAAACTCAACTTGTTACGTCTTGGGTGACCGTAAGATTTAGGTCCGTCATCTCGGTGGTAGTGGAGCTGGCCCACATTCGTT"
                .as_bytes()
                .to_owned();
            let base_qualities = vec![40; pattern.len()];

            let mut stack = MinMaxHeap::new();
            let mut tree = Tree::new();
            let intervals = k_mismatch_search(
                &pattern,
                &base_qualities,
                &parameters,
                &fmd_index,
                &mut stack,
                &mut tree,
                &difference_model,
                &mismatch_bound,
            );
            assert_eq!(intervals.len(), 1);
        }
    }
}
