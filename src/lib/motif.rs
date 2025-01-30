use crate::{iupac::IupacBase, modtype::ModType};
use anyhow::{Ok, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComplementMotif {
    pub sequence: Vec<IupacBase>,
    pub mod_type: ModType,
    pub position: u8,
}

pub trait MotifLike {
    fn sequence_string(&self) -> String;
    fn regex(&self) -> Result<String>;
    fn as_string(&self) -> String;
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Motif {
    pub sequence: Vec<IupacBase>,
    pub mod_type: ModType,
    pub position: u8,
}

impl Motif {
    pub fn new(sequence: &str, mod_type: &str, position: u8) -> Result<Self, anyhow::Error> {
        let parsed_sequence = sequence
            .chars()
            .map(|c| {
                IupacBase::from_char(c).map_err(|_| {
                    anyhow::anyhow!("Error parsing base: {}, in sequence: {}. ", c, sequence)
                })
            })
            .collect::<Result<Vec<IupacBase>, anyhow::Error>>()?;
        let mod_type = mod_type.parse::<ModType>()?;

        let base_at_position = parsed_sequence.get(position as usize).ok_or_else(|| {
            anyhow::anyhow!(
                "Position {} is out of bounds for sequence: {}",
                position,
                sequence
            )
        })?;
        if base_at_position != mod_type.get_iupac_base() {
            anyhow::bail!(
                "Base at position {} ({}) does not match mod type: {}",
                position,
                base_at_position.to_string(),
                mod_type.to_string()
            );
        }
        Ok(Self {
            sequence: parsed_sequence,
            mod_type,
            position,
        })
    }

    pub fn reverse_complement(&self) -> Result<ComplementMotif, anyhow::Error> {
        let reversed_sequence: String = self
            .sequence
            .iter()
            .rev()
            .map(|base| base.complement().to_string())
            .collect::<String>();
        let position = self.sequence.len() as u8 - self.position - 1;

        let mod_type_str: &str = &self.mod_type.to_string();
        ComplementMotif::new(&reversed_sequence, mod_type_str, position)
    }

    pub fn reverse_complement_sequence(&self) -> String {
        let reversed_sequence: String = self
            .sequence
            .iter()
            .rev()
            .map(|base| base.complement().to_string())
            .collect::<String>();
        reversed_sequence
    }

    // Function that grows a motif at a given position with a set of bases which can be outside the motif sequence e.g. GATC to CNGATC

    pub fn add_base_upstream(&mut self, base: IupacBase) -> Result<&Self> {
        self.sequence.insert(0, base);
        self.position += 1;
        Ok(self)
    }

    pub fn add_base_downstream(&mut self, base: IupacBase) -> Result<&Self> {
        self.sequence.push(base);
        Ok(self)
    }

    pub fn grow_motif_at(&mut self, position: i8, bases: Vec<&str>) -> Result<&Self> {
        let new_iupac = IupacBase::from_set(bases)?;
        match position {
            // if ngeative grow upstream
            p if p < 0 => {
                let p = p.abs();
                for _ in 0..(p - 1) {
                    self.add_base_upstream(IupacBase::N)?;
                }
                self.add_base_upstream(new_iupac)?;
                Ok(self)
            }
            // if positive grow downstream
            p if p >= self.sequence.len() as i8 => {
                for _ in 0..(p - self.sequence.len() as i8) {
                    self.add_base_downstream(IupacBase::N)?;
                }
                self.add_base_downstream(new_iupac)?;
                Ok(self)
            }
            // if position is within the sequence
            _ => {
                let old_iupac = self.sequence[position as usize];
                let updated_iupac = new_iupac.join_other(old_iupac)?;

                self.sequence[position as usize] = updated_iupac;
                Ok(self)
            }
        }
    }
    pub fn as_pretty_string(&self) -> String {
        // replace iupac string modification type at the mod position
        let mut pretty_sequence = self.sequence_string();
        pretty_sequence.replace_range(
            self.position as usize..(self.position + 1) as usize,
            &self.mod_type.to_string(),
        );

        // return the pretty string
        pretty_sequence
    }
}
impl MotifLike for Motif {
    fn sequence_string(&self) -> String {
        self.sequence.iter().map(|b| b.to_string()).collect()
    }

    fn regex(&self) -> Result<String> {
        let mut regex = String::new();
        for base in self.sequence.iter() {
            regex.push_str(&base.to_regex());
        }
        Ok(regex)
    }

    fn as_string(&self) -> String {
        let sequence: String = self.sequence.iter().map(|b| b.to_string()).collect();
        format!(
            "{}_{}_{}",
            sequence,
            self.mod_type.to_string(),
            self.position
        )
    }


}

impl ComplementMotif {
    pub fn new(sequence: &str, mod_type: &str, position: u8) -> Result<Self, anyhow::Error> {
        let parsed_sequence = sequence
            .chars()
            .map(|c| {
                IupacBase::from_char(c).map_err(|_| {
                    anyhow::anyhow!("Error parsing base: {}, in sequence: {}. ", c, sequence)
                })
            })
            .collect::<Result<Vec<IupacBase>, anyhow::Error>>()?;
        let mod_type = mod_type.parse::<ModType>()?;

        let base_at_position = parsed_sequence.get(position as usize).ok_or_else(|| {
            anyhow::anyhow!(
                "Position {} is out of bounds for sequence: {}",
                position,
                sequence
            )
        })?;
        if base_at_position != &mod_type.get_iupac_base().complement() {
            anyhow::bail!(
                "The complement base at position {} ({}) does not match mod type: {}",
                position,
                base_at_position.complement().to_string(),
                mod_type.to_string()
            );
        }
        Ok(Self {
            sequence: parsed_sequence,
            mod_type,
            position,
        })
    }
}

impl MotifLike for ComplementMotif {
    fn sequence_string(&self) -> String {
        self.sequence.iter().map(|b| b.to_string()).collect()
    }

    fn regex(&self) -> Result<String> {
        let mut regex = String::new();
        for base in self.sequence.iter() {
            regex.push_str(&base.to_regex());
        }
        Ok(regex)
    }

    fn as_string(&self) -> String {
        let sequence: String = self.sequence.iter().map(|b| b.to_string()).collect();
        format!(
            "{}_{}_{}",
            sequence,
            self.mod_type.to_string(),
            self.position
        )
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MotifPair {
    pub forward: Motif,
    pub reverse: Motif,
    pub is_palindromic: bool,
}

impl MotifPair {
    pub fn new(forward: Motif, reverse: Motif) -> Result<Self, anyhow::Error> {
        if forward.sequence.len() != reverse.sequence.len() {
            anyhow::bail!("Forward and reverse motifs must have the same length");
        }
        // Check that the forward and reverse motifs are reverse complements of each other
        for (f_base, r_base) in forward.sequence.iter().zip(reverse.sequence.iter().rev()) {
            if f_base != &r_base.complement() {
                anyhow::bail!(
                    "Forward and reverse motifs are not reverse complements of each other"
                );
            }
        }
        let is_palindromic = forward == reverse;
        Ok(Self {
            forward,
            reverse,
            is_palindromic,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        assert_eq!(
            motif.sequence,
            vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::T]
        );
        assert_eq!(motif.mod_type, ModType::SixMA);
        assert_eq!(motif.position, 0);

        let motif = Motif::new("GATC", "6mA", 1);
        assert!(motif.is_ok());
        let motif = Motif::new("GATC", "5mC", 3);
        assert!(motif.is_ok());
        let motif = Motif::new("GATC", "4mC", 3);
        assert!(motif.is_ok());

        // Wrong base at position
        let motif = Motif::new("ACGT", "6mA", 1);
        assert!(motif.is_err());
        // Wrong base at position
        let motif = Motif::new("ACGT", "m", 0);
        assert!(motif.is_err());
        // Wrong base at position
        let motif = Motif::new("ACGT", "5mC", 0);
        assert!(motif.is_err());
        // Position out of bounds
        let motif = Motif::new("ACGT", "a", 4);
        assert!(motif.is_err());
        // Position out of bounds
        let motif = Motif::new("ACGT", "a", 1);
        assert!(motif.is_err());
        // Invalid sequence
        let motif = Motif::new("LKOIUASUJDNG", "5mC", 0);
        assert!(motif.is_err());
        // Invalid sequence
        let motif = Motif::new("GATC.", "5mC", 0);
        assert!(motif.is_err());
        // Invalid modtype
        let motif = Motif::new("GATC", "6cA", 1);
        assert!(motif.is_err());
    }

    #[test]
    fn test_motif_grow_motif_at() {
        let motif_default = Motif::new("ACGT", "6mA", 0).unwrap();
        let mut motif = motif_default.clone();
        motif.grow_motif_at(0, vec!["A", "C"]).unwrap();
        assert_eq!(
            motif.sequence,
            vec![IupacBase::M, IupacBase::C, IupacBase::G, IupacBase::T]
        );
        let mut motif = motif_default.clone();
        motif.grow_motif_at(1, vec!["A", "C"]).unwrap();
        assert_eq!(
            motif.sequence,
            vec![IupacBase::A, IupacBase::M, IupacBase::G, IupacBase::T]
        );
        let mut motif = motif_default.clone();
        motif.grow_motif_at(2, vec!["A", "C"]).unwrap();
        assert_eq!(
            motif.sequence,
            vec![IupacBase::A, IupacBase::C, IupacBase::V, IupacBase::T]
        );
        let mut motif = motif_default.clone();
        motif.grow_motif_at(3, vec!["A", "C"]).unwrap();
        assert_eq!(
            motif.sequence,
            vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::H]
        );
    }

    #[test]
    fn test_reverse_complement() {
        let motif = Motif::new("AACT", "6mA", 0).unwrap();
        let revcomp = motif.reverse_complement().unwrap();
        assert_eq!(
            revcomp.sequence,
            vec![IupacBase::A, IupacBase::G, IupacBase::T, IupacBase::T]
        );
        assert_eq!(revcomp.mod_type, ModType::SixMA);
        assert_eq!(revcomp.position, 3);

        let motif = Motif::new("ACGT", "5mC", 1).unwrap();
        let revcomp = motif.reverse_complement().unwrap();
        assert_eq!(
            revcomp.sequence,
            vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::T]
        );
        assert_eq!(revcomp.mod_type, ModType::FiveMC);
        assert_eq!(revcomp.position, 2);

        let motif = Motif::new("ATAT", "6mA", 2).unwrap();
        let revcomp = motif.reverse_complement().unwrap();
        assert_eq!(
            revcomp.sequence,
            vec![IupacBase::A, IupacBase::T, IupacBase::A, IupacBase::T]
        );
        assert_eq!(revcomp.mod_type, ModType::SixMA);
        assert_eq!(revcomp.position, 1);

        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let revcomp = motif.reverse_complement().unwrap();
        assert_eq!(
            revcomp.sequence,
            vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::T]
        );
        assert_eq!(revcomp.mod_type, ModType::SixMA);
        assert_eq!(revcomp.position, 3);

        let motif = Motif::new("CAGCTG", "4mC", 3).unwrap();
        let revcomp = motif.reverse_complement().unwrap();
        assert_eq!(
            revcomp.sequence,
            vec![
                IupacBase::C,
                IupacBase::A,
                IupacBase::G,
                IupacBase::C,
                IupacBase::T,
                IupacBase::G
            ]
        );
        assert_eq!(revcomp.mod_type, ModType::FourMC);
        assert_eq!(revcomp.position, 2);
    }

    #[test]
    fn test_motif_pair_new() {
        let forward = Motif::new("CCANNGT", "6mA", 2).unwrap();
        let reverse = Motif::new("ACNNTGG", "6mA", 0).unwrap();
        let pair = MotifPair::new(forward, reverse);
        assert!(pair.is_ok());
        // Different length error
        let forward = Motif::new("ACGT", "6mA", 0).unwrap();
        let reverse = Motif::new("ACGTT", "6mA", 0).unwrap();
        let pair = MotifPair::new(forward, reverse);
        assert!(pair.is_err());
        // Not reverse complement error
        let forward = Motif::new("ACGT", "6mA", 0).unwrap();
        let reverse = Motif::new("AAAT", "6mA", 2).unwrap();
        let pair = MotifPair::new(forward, reverse);
        assert!(pair.is_err());
    }
}
