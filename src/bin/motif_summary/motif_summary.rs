use log::{info, debug, warn};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use rustymotif_utils::motif;
use rustymotif_utils::pileup;
use rustymotif_utils::modtype::ModType;
use rustymotif_utils::strand::Strand;
use std::time::Instant;
use std::fs::File;
use csv::{WriterBuilder};
use anyhow::Result;
use crate::data::AnnotatedPileupRecord;
use crate::io::{self, load_bin_motifs};
use crate::io::load_contig_bin;
use crate::io::AnnotatedRecordWriter;
use crate::cli;
use crate::fasta_reader;
use crate::sequence;
use crate::data;

pub fn motif_summary(args: &cli::Cli) -> Result<()> {
    let global_timer = Instant::now();

    info!("Starting motif methylation state processing");
    info!("Loading reference file: {}", args.reference);
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let bed_file = File::open(&args.bed)
        .map_err(|e| anyhow::anyhow!("Could not open bed file: {} ({})", args.bed, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(bed_file, 0);
   
    let contig_bin_path = Path::new(&args.contig_bin);
    info!("Loading contig bin from: {}", args.contig_bin);
    let bin_contig_map = load_contig_bin(contig_bin_path)?;

    let bin_motifs_path = Path::new(&args.bin_motifs);
    info!("Loading bin motifs from: {}", args.bin_motifs);
    let mut bin_motifs = load_bin_motifs(bin_motifs_path)
        .map_err(|e| anyhow::anyhow!("Error loading bin motifs: {}", e))?;
    debug!("Loaded {} bin motifs", bin_motifs.values().len());


    // Create a map to hold motifs for each contig
    let mut contig_motif_map: HashMap<String, Vec<motif::Motif>> = HashMap::new();
    for (bin, motifs) in bin_motifs.iter_mut() {
        if let Some(contigs) = bin_contig_map.get(bin) {
            for contig in contigs {
                contig_motif_map.entry(contig.clone())
                    .or_insert_with(Vec::new)
                    .extend(motifs.iter().cloned());
            }
        } else {
            warn!("No contigs found for bin: {}", bin);
        }
    }

    let outdir = Path::new(&args.out);
    let output: &Path = &outdir.join(format!("motif_annotated_pileup.bed"));
    let mut result_writer = AnnotatedRecordWriter::new(output)?;
    result_writer.write_header()?;

    let summary_output: &Path = &outdir.join(format!("motif_annotated_pileup_summary.tsv"));
    let mut summary_writer = io::AnnotatedRecordSummaryWriter::new(summary_output)?;
    summary_writer.write_header()?;
    
    info!("Processing pileup file: {}", args.bed);
    loop {
        info!("Processing a batch");
        let timer = Instant::now();
        let chunks = pileup_reader.load_n_chunks(1);
        match chunks {
            Some(chunks) => {
                info!("Loaded batch {:?}", timer.elapsed());
                let mut builder = data::GenomeWorkSpaceBuilder::new();
                for chunk in chunks {
                    let contig_id = &chunk.reference;
                    info!("Processing contig: {}", contig_id);
                    debug!("Adding contig to workspace");
                    let contig = sequence::Contig::new(
                        contig_id,
                        reference.get(contig_id).unwrap().as_str(),
                    );
                    builder.add_contig(contig);
                    debug!("Adding records to contig");
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                
                
                for ( contig_id, contig) in genome_work_space.contigs.iter() {
                    debug!("Contig ID: {}", contig_id);
                    let motifs = contig_motif_map.get(contig_id)
                        .cloned()
                        .unwrap_or_else(|| vec![]);
                    let motifs_set: HashSet<motif::Motif> = HashSet::from_iter(motifs.iter().cloned());

                    // Iterate over pileup records in the contig
                    for mod_type in ModType::iter() {
                        debug!("Processing mod_type: {}", mod_type);
                        let records = contig.get_by_mod_type(mod_type.clone());
                        if records.is_empty() {
                            continue;
                        }

                        let mut annotated_records = records.iter()
                            .map(|record| {
                                AnnotatedPileupRecord::new(
                                    record.clone().clone(),
                                    false, // in_motif will be set later
                                    vec![], // Clone motifs for each record
                                )
                            })
                            .collect::<Vec<AnnotatedPileupRecord>>();
                        
        
                        // Add motif to records where motif occours
                        for motif in motifs_set.iter() {
                            if motif.mod_type != mod_type {
                                continue;
                            }
                            debug!("Processing motif: {}", motif.as_pretty_string());
                        
                            // Get index of where motif occours
                            let motif_index = contig.find_motif_indices(motif)?;
                            let motif_index_set = match motif_index {
                                None => None,
                                Some(indices) => Some(indices.iter().cloned().collect::<HashSet<usize>>()),
                            };

                            let motif_index_reverse = contig.find_complement_motif_indices(motif)?;
                            let motif_index_set_reverse = match motif_index_reverse {
                                None => None,
                                Some(indices) => Some(indices.iter().cloned().collect::<HashSet<usize>>()),
                            };


                            for record in annotated_records.iter_mut() {
                                match record.pileup_record.strand {
                                    Strand::Positive => {
                                        match &motif_index_set {
                                            None => continue,
                                            Some(set) => {
                                                if set.contains(&record.pileup_record.position) {
                                                    record.in_motif = true;
                                                    record.motifs.push(motif.clone());
                                                }
                                            }
                                        }
                                    },
                                    Strand::Negative => {
                                        match &motif_index_set_reverse {
                                            None => continue,
                                            Some(set) => {
                                                if set.contains(&record.pileup_record.position) {
                                                    record.in_motif = true;
                                                    record.motifs.push(motif.clone());
                                                }
                                            }
                                        }
                                    }
                                }
                            }
    
                        }
                        // Write annotated records to file
                        result_writer.write_records_iter(&annotated_records)?;
                        let mut summary = data::AnnotatedPileupSummary::new();
                        summary.summarise_records(
                            &annotated_records,
                            args.summary_min_cov,
                            args.summary_min_methylation,
                        )?;
                        summary_writer.write_summary(
                            summary,
                        )?;
                    }   
                }
            }
            None => {
                info!("Contig did not contain any records");
            }
        }

        info!("Finished batch in {:?}", timer.elapsed());

        if pileup_reader.eof_reached {
            break;
        }
    }
    info!("Finished processing in {:?}", global_timer.elapsed());
    Ok(())
}