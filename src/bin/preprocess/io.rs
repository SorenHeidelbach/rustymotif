
use std::path::Path;
use std::fs::File;
use std::fs;
use anyhow::Result;
use csv::{WriterBuilder, Writer};
use rustymotif_utils::{motif::Motif, pileup::PileupRecord};
use log::{info, debug, warn};
use rust_htslib::bgzf::Writer as BgzfWriter;
use std::io::Write;
use std::io;
use std::path::PathBuf;


pub trait MethylBedWriter {
    fn write_header(&mut self) -> Result<()>;
    fn write_record(&mut self, record: &PileupRecord) -> Result<()>;
    fn flush(&mut self) -> Result<()>;
}

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
impl MethylBedWriter for RecordWriter {
    fn write_header(&mut self) -> Result<()> {
        let methyl_bed_fields = rustymotif_utils::pileup::MethylBedFields::all_as_string();
        self.writer.write_record(&methyl_bed_fields)?;
        self.writer.flush()?;
        Ok(())
    }

    fn write_record(&mut self, record: &PileupRecord) -> Result<()> {
        let fields = record.to_bed_fields()?;
        self.writer.write_record(&fields)?;
        Ok(())
    }
    
    fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        self.file.sync_all()?;
        Ok(())
    }

}


pub struct BgzippedRecordWriter {
    writer: BgzfWriter,
}

impl BgzippedRecordWriter {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let writer = BgzfWriter::from_path(path)?;
        Ok(Self { writer })
    }
    pub fn write_records_iter<'a, I>(&mut self, records: I) -> Result<()>
    where
        I: IntoIterator<Item = &'a PileupRecord>,
    {
        for record in records {
            self.write_record(record)?;
        }
        self.flush()?;
        Ok(())
    }
}
impl MethylBedWriter for BgzippedRecordWriter {
    fn write_header(&mut self) -> Result<()> {
        let fields = rustymotif_utils::pileup::MethylBedFields::all_as_string();
        let header = fields
            .iter()
            .map(|field| field.to_string())
            .collect::<Vec<String>>()
            .join("\t");

        self.writer.write_all(header.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    fn write_record(&mut self, record: &PileupRecord) -> Result<()> {
        let line = record.to_bed_fields()?.join("\t");
        self.writer.write_all(line.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }


}


pub fn create_methyl_bed_writer<P: AsRef<Path>>(path: P, bgzip: bool) -> Result<Box<dyn MethylBedWriter>> {
    let path = path.as_ref();

    if bgzip {
        debug!("Output will be written in bgzipped format to: {}", path.display());

        // Check and possibly adjust extension
        let final_path: PathBuf = if path.extension().map_or(true, |ext| ext != "gz") {
            warn!("Output file does not have .gz extension, appending .gz automatically.");
            let new_path = path.with_extension("gz");
            info!("Adjusted output file path: {}", new_path.display());
            new_path
        } else {
            path.to_path_buf()
        };

        Ok(Box::new(BgzippedRecordWriter::from_path(&final_path)?))
    } else {
        Ok(Box::new(RecordWriter::new(path)?))
    }
}

pub fn write_records_iter<'a, I>(
    writer: &mut dyn MethylBedWriter,
    records: I,
) -> Result<()>
where
    I: IntoIterator<Item = &'a PileupRecord>,
{
    for record in records {
        writer.write_record(record)?;
    }
    writer.flush()
}