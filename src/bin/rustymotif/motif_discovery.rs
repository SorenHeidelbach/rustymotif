use crate::{cli, data, fasta_reader, search, sequence};
use anyhow::{bail, Result};
use bio::bio_types::annot::contig;
use log::{debug, info};
use std::{fs::File, path::Path, time::Instant};
use rustymotif_utils::{motif, motif::MotifLike, pileup, strand::Strand};
use csv::{WriterBuilder};

pub fn rustymotif(args: &cli::Cli) -> Result<()> {
    let global_timer = Instant::now();
    info!("Loading reference file: {}", args.reference);
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
                        &sequence::MethylationThresholds::new(args.middle_threshold as f32, args.high_threshold as f32),
                    );
                    builder.add_contig(contig);
                    debug!("Adding records to contig");
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                for (refenrece_id, contig) in genome_work_space.contigs.into_iter() {
                    println!("Processing contig: {}", refenrece_id);
                    let identified_motifs = search::motif_search(
                        contig,
                        args.window_size,
                        0.05,
                        0.2,
                        0.25,
                        15,
                        4
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

