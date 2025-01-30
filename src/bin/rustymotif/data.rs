use crate::pileup::PileupChunk;
use crate::sequence::Contig;
use ahash::{HashMap, HashMapExt};

pub struct GenomeWorkSpaceBuilder {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkSpaceBuilder {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: Contig) {
        self.contigs.insert(contig.reference.clone(), contig);
    }

    pub fn push_records(&mut self, records: PileupChunk) {
        let contig = self
            .contigs
            .get_mut(records.reference.as_str())
            .expect("Could not find contig");
        contig.add_records(records);
    }

    pub fn build(self) -> GenomeWorkspace {
        GenomeWorkspace {
            contigs: self.contigs,
        }
    }
}

pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}
