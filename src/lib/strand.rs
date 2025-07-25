use anyhow::{bail, Result};
use std::{fmt::Display, str::FromStr};
use serde::Serialize;

/// Represents the DNA strand of reference.
#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy, Serialize)]
pub enum Strand {
    Positive,
    Negative,
}

impl Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}


impl Strand {
    pub fn to_string(&self) -> String {
        match self {
            Strand::Positive => "+".to_string(),
            Strand::Negative => "-".to_string(),
        }
    }
}

impl FromStr for Strand {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Positive),
            "-" => Ok(Strand::Negative),
            _ => bail!("Could not parse '{}' to Strand", s),
        }
    }
}
