use crate::{cli, data, fasta_reader, search, sequence, score_motifs};
use anyhow::{bail, Result};
use bio::bio_types::annot::contig;
use log::{debug, info};
use std::{fs::File, path::Path, time::Instant};
use utils::{motif, motif::MotifLike, pileup, strand::Strand};
use csv::{WriterBuilder};

pub fn rustymotif(args: &cli::Cli) -> Result<()> {
    let global_timer = Instant::now();
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let pileup_file = File::open(&args.pileup)
        .map_err(|e| anyhow::anyhow!("Could not open pileup file: {} ({})", args.pileup, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(pileup_file, args.min_cov);

    let outdir = Path::new(&args.out);
    let out_motif_path = outdir.join("motifs.tsv");
    let out_motif_file = File::create(out_motif_path)
        .map_err(|e| anyhow::anyhow!("Could not create output file: {} ({})", args.out, e))?;
    let mut motif_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(out_motif_file);

    info!("Processing pileup file: {}", args.pileup);
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

                let mut motifs = Vec::new();
                motifs.push(
                    motif::Motif::new(
                        "RGATCY",
                        "4mC",
                        4,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "GATC",
                        "6mA",
                        1,
                    )?
                );
                
                motifs.push(
                    motif::Motif::new(
                        "A",
                        "6mA",
                        0,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "GA",
                        "6mA",
                        1,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "CGA",
                        "6mA",
                        2,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "TCGA",
                        "6mA",
                        3,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "CTCGA",
                        "6mA",
                        4,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "TCGAG",
                        "6mA",
                        3,
                    )?
                );
                motifs.push(
                    motif::Motif::new(
                        "CTCGAG",
                        "6mA",
                        4,
                    )?
                );
                // Write motif_beta_true, motif_mixture, motif_log_likelihood, motif_log_likelihood_false to a file
                let output = outdir.join(format!("mixed_model_of_motifs.tsv"));
                let output_string = output.to_str()
                    .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;
                
                let result_file = File::create(output_string)?;
                let mut result_writer = WriterBuilder::new()
                    .delimiter(b'\t')
                    .has_headers(true)
                    .from_writer(result_file);
                result_writer.write_record(&[
                    "contig_id",
                    "motif",
                    "beta_methylated_shape1",
                    "beta_methylated_shape2",
                    "beta_unmethylated_shape1",
                    "beta_unmethylated_shape2",
                    "mixture",
                    "log_likelihood_methylated",
                    "log_likelihood_unmethylated",
                ])?;
                for ( contig_id, contig) in genome_work_space.contigs.iter() {
                    for motif in motifs.iter() {
                        debug!("Searching for motif: {}", motif.as_pretty_string());
                        debug!("Contig ID: {}", contig_id);
                        // Make output String
                        let output = outdir.join(format!("{}_motif_{}.tsv", contig_id, motif.as_pretty_string()));
                        let output_string = output.to_str()
                            .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;
                        
                        let motif_mixture_model = match score_motifs::score_motif(
                            motif.clone(),
                            contig.clone(),
                            Some(output_string.to_string()),
                        )? {
                            Some(motif_mixture_model) => {
                                debug!("Motif found in contig: {}", contig_id);
                                motif_mixture_model
                            }
                            None => {
                                debug!("Motif not found in contig: {}", contig_id);
                                continue;
                            }
                        };
                        debug!("Motif search complete for contig: {}", contig_id);

                        let (log_lik_true, log_lik_false) = motif_mixture_model.log_likelihoods(
                            &contig.get_by_mod_type(motif.mod_type.clone()),
                        );
                        result_writer.write_record(&[
                            contig_id,
                            &motif.as_pretty_string(),
                            &motif_mixture_model.true_model.alpha.to_string(),
                            &motif_mixture_model.true_model.beta.to_string(),
                            &motif_mixture_model.false_model.alpha.to_string(),
                            &motif_mixture_model.false_model.beta.to_string(),
                            &motif_mixture_model.pi.to_string(),
                            &log_lik_true.to_string(),
                            &log_lik_false.to_string(),
                        ])?;



                        
                    }
                }
                result_writer.flush()?;

                for (refenrece_id, contig) in genome_work_space.contigs.into_iter() {
                    println!("Processing contig: {}", refenrece_id);
                    let identified_motifs = search::motif_search(
                        contig,
                        10,
                        50,
                        0.01,
                        Some("intermediate_motifs.tsv"),
                    )?;
                    for motif in identified_motifs {
                        motif_writer.serialize(motif)?;
                        motif_writer.flush()?;
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

