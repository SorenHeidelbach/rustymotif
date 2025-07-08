
use std::path::Path;
use std::fs::File;
use std::fs;
use anyhow::Result;
use csv::{WriterBuilder, Writer};
use rustymotif_utils::{motif::Motif, pileup::PileupRecord};
use log::{info, debug, warn};




pub struct RecordWriter {
    file: File,
    writer: Writer<File>,
}
impl RecordWriter {
    pub fn new(file_path: &Path) -> Result<Self> {
        let file = File::create(file_path)?;
        let writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(file.try_clone()?);
        Ok(Self { file, writer })
    }
    pub fn write_header(&mut self) -> Result<()> {
        self.writer.write_record(&[
            "chrom",               // BED column 1
            "start position",      // BED column 2
            "end position",        // BED column 3
            "modified base code",  // BED column 4 (name field)
            "score",               // BED column 5
            "strand",              // BED column 6
            "start position",      // compatibility
            "end position",        // compatibility
            "color",               // BED RGB
            "Nvalid_cov",
            "fraction modified",
            "Nmod",
            "Ncanonical",
            "Nother_mod",
            "Ndelete",
            "Nfail",
            "Ndiff",
            "Nnocall",
        ])?;
        self.writer.flush()?;
        Ok(())
    }

    pub fn write_record(&mut self, record: &PileupRecord) -> Result<()> {
        let pileup = &record;

        let chrom = &pileup.reference;
        let start = pileup.position;
        let end = start + 1; // BED uses exclusive end
        let mod_code = pileup.mod_type.to_pileup_code().to_owned(); // assuming this gives a single-letter code
        let score = pileup.n_valid_cov;
        let strand = pileup.strand.to_string();

        let fraction_modified = if pileup.n_valid_cov > 0 {
            format!("{:.3}", pileup.n_mod as f32 / pileup.n_valid_cov as f32)
        } else {
            "0.000".to_string()
        };

        self.writer.write_record(&[
            chrom,
            &start.to_string(),
            &end.to_string(),
            &mod_code,
            &score.to_string(),
            &strand,
            &start.to_string(),     // compatibility
            &end.to_string(),       // compatibility
            "255,0,0",              // BED RGB field
            &pileup.n_valid_cov.to_string(),
            &fraction_modified,
            &pileup.n_mod.to_string(),
            &pileup.n_canonical.to_string(),
            "NA",
            "NA",
            "NA",
            &pileup.n_diff.to_string(),
            "NA",
        ])?;
        Ok(())
    }
    
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        self.file.sync_all()?;
        Ok(())
    }

    pub fn write_records_iter<'a, I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = &'a PileupRecord>,
    {
        for record in records.into_iter() {
            self.write_record(record)?;
        }
        self.flush()?;
        Ok(())
    }
}
