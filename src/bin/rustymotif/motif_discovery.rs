use crate::{cli, data, fasta_reader, search, sequence};
use anyhow::{bail, Result};
use log::{debug, info};
use std::{fs::File, path::Path, time::Instant};
use utils::{motif, motif::MotifLike, pileup, strand::Strand};

pub fn rustymotif(args: &cli::Cli) -> Result<(), anyhow::Error> {
    let global_timer = Instant::now();
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let pileup_file = File::open(&args.pileup)
        .map_err(|e| anyhow::anyhow!("Could not open pileup file: {} ({})", args.pileup, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(pileup_file, args.min_cov);

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
                    builder.add_contig(reference.get(contig_id).unwrap().clone());
                    debug!("Adding records to contig");
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                for (refenrece_id, contig) in genome_work_space.contigs.into_iter() {
                    println!("Processing contig: {}", refenrece_id);
                    search::motif_search(
                        contig, 
                        args.window_size, 
                        args.min_kl_divergence,
                        0.25,
                        10
                    )
                        .map_err(|e| anyhow::anyhow!("Error searching for motifs: {}", e))?;
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
