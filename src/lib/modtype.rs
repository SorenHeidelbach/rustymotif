use crate::iupac::IupacBase;
use anyhow::{bail, Result};
use std::{fmt, str::FromStr};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum ModType {
    SixMA,
    FiveMC,
    FourMC,
}

impl ModType {
    pub fn to_pileup_code(&self) -> &'static str {
        match self {
            ModType::SixMA => "a",
            ModType::FiveMC => "m",
            ModType::FourMC => "21839",
        }
    }

    pub fn get_iupac_base(&self) -> &IupacBase {
        match self {
            ModType::SixMA => &IupacBase::A,
            ModType::FiveMC => &IupacBase::C,
            ModType::FourMC => &IupacBase::C,
        }
    }

    pub fn to_string(&self) -> &'static str {
        match self {
            ModType::SixMA => "6mA",
            ModType::FiveMC => "5mC",
            ModType::FourMC => "4mC",
        }
    }
}

impl fmt::Display for ModType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ModType::FiveMC => write!(f, "5mC, (m)"),
            ModType::FourMC => write!(f, "4mC, (21839)"),
            ModType::SixMA => write!(f, "6mA, (a)"),
        }
    }
}

impl FromStr for ModType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "5mC" | "m" => Ok(ModType::FiveMC),
            "4mC" | "21839" => Ok(ModType::FourMC),
            "6mA" | "a" => Ok(ModType::SixMA),
            _ => bail!("Invalid ModType: {}", s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_pileup_code() {
        for mt in vec![ModType::SixMA, ModType::FiveMC, ModType::FourMC] {
            let pileup_code = mt.to_pileup_code();
            let expected = match mt {
                ModType::SixMA => "a",
                ModType::FiveMC => "m",
                ModType::FourMC => "21839",
            };
            assert_eq!(pileup_code, expected)
        }
    }

    #[test]
    fn test_display() {
        for mt in vec![ModType::SixMA, ModType::FiveMC, ModType::FourMC] {
            let display = format!("{}", mt);
            let expected = match mt {
                ModType::SixMA => "6mA, (a)",
                ModType::FiveMC => "5mC, (m)",
                ModType::FourMC => "4mC, (21839)",
            };
            assert_eq!(display, expected)
        }
    }

    #[test]
    fn test_from_str() {
        for (s, expected) in vec![
            ("5mC", ModType::FiveMC),
            ("m", ModType::FiveMC),
            ("4mC", ModType::FourMC),
            ("21839", ModType::FourMC),
            ("6mA", ModType::SixMA),
            ("a", ModType::SixMA),
        ] {
            let modtype = s.parse::<ModType>().unwrap();
            assert_eq!(modtype, expected)
        }
    }
    #[test]
    fn test_from_str_error() {
        for s in vec!["5mc", "4mc", "6ma", "b"] {
            let modtype = s.parse::<ModType>();
            assert!(modtype.is_err())
        }
    }
}
