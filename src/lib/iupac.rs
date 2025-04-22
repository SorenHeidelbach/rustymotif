use anyhow::bail;
use anyhow::Result;
use std::fmt::Display;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, EnumIter)]
pub enum IupacBase {
    A,
    C,
    G,
    T,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    N,
}

impl Display for IupacBase {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

impl IupacBase {
    pub fn to_regex(&self) -> &'static str {
        match self {
            IupacBase::A => "A",
            IupacBase::C => "C",
            IupacBase::G => "G",
            IupacBase::T => "T",
            IupacBase::R => "[AG]",
            IupacBase::Y => "[CT]",
            IupacBase::S => "[GC]",
            IupacBase::W => "[AT]",
            IupacBase::K => "[GT]",
            IupacBase::M => "[AC]",
            IupacBase::B => "[CGT]",
            IupacBase::D => "[AGT]",
            IupacBase::H => "[ACT]",
            IupacBase::V => "[ACG]",
            IupacBase::N => ".",
        }
    }

    pub fn to_set(&self) -> Vec<&'static str> {
        match self {
            IupacBase::A => vec!["A"],
            IupacBase::C => vec!["C"],
            IupacBase::G => vec!["G"],
            IupacBase::T => vec!["T"],
            IupacBase::R => vec!["A", "G"],
            IupacBase::Y => vec!["C", "T"],
            IupacBase::S => vec!["G", "C"],
            IupacBase::W => vec!["A", "T"],
            IupacBase::K => vec!["G", "T"],
            IupacBase::M => vec!["A", "C"],
            IupacBase::B => vec!["C", "G", "T"],
            IupacBase::D => vec!["A", "G", "T"],
            IupacBase::H => vec!["A", "C", "T"],
            IupacBase::V => vec!["A", "C", "G"],
            IupacBase::N => vec!["A", "C", "G", "T"],
        }
    }

    pub fn join_other(&self, other: IupacBase) -> Result<IupacBase> {
        let mut set = self.to_set();
        set.extend(other.to_set());
        IupacBase::from_set(set)
    }

    pub fn from_set(set: Vec<&str>) -> Result<IupacBase> {
        // keep only unique elements
        let set = set
            .into_iter()
            .collect::<std::collections::HashSet<&str>>()
            .into_iter()
            .collect::<Vec<&str>>();
        match set.len() {
            1 => match set[0] {
                "A" => Ok(IupacBase::A),
                "C" => Ok(IupacBase::C),
                "G" => Ok(IupacBase::G),
                "T" => Ok(IupacBase::T),
                _ => bail!("Invalid IUPAC base: {}", set[0]),
            },
            2 => {
                let mut set = set.iter().map(|&x| x).collect::<Vec<&str>>();
                set.sort();
                match set.as_slice() {
                    ["A", "G"] => Ok(IupacBase::R),
                    ["C", "T"] => Ok(IupacBase::Y),
                    ["C", "G"] => Ok(IupacBase::S),
                    ["A", "T"] => Ok(IupacBase::W),
                    ["G", "T"] => Ok(IupacBase::K),
                    ["A", "C"] => Ok(IupacBase::M),
                    _ => bail!("Invalid IUPAC base: {:?}", set),
                }
            }
            3 => {
                let mut set = set.iter().map(|&x| x).collect::<Vec<&str>>();
                set.sort();
                match set.as_slice() {
                    ["C", "G", "T"] => Ok(IupacBase::B),
                    ["A", "G", "T"] => Ok(IupacBase::D),
                    ["A", "C", "T"] => Ok(IupacBase::H),
                    ["A", "C", "G"] => Ok(IupacBase::V),
                    _ => bail!("Invalid IUPAC base: {:?}", set),
                }
            }
            4 => {
                let mut set = set.iter().map(|&x| x).collect::<Vec<&str>>();
                set.sort();
                match set.as_slice() {
                    ["A", "C", "G", "T"] => Ok(IupacBase::N),
                    _ => bail!("Invalid IUPAC base: {:?}", set),
                }
            }
            _ => bail!("Invalid IUPAC base: {:?}", set),
        }
    }

    pub fn to_onehot(&self) -> [f64; 4] {
        match self {
            IupacBase::A => [1.0, 0.0, 0.0, 0.0],
            IupacBase::C => [0.0, 1.0, 0.0, 0.0],
            IupacBase::G => [0.0, 0.0, 1.0, 0.0],
            IupacBase::T => [0.0, 0.0, 0.0, 1.0],
            IupacBase::R => [0.5, 0.0, 0.5, 0.0],
            IupacBase::Y => [0.0, 0.5, 0.0, 0.5],
            IupacBase::S => [0.0, 0.5, 0.5, 0.0],
            IupacBase::W => [0.5, 0.0, 0.0, 0.5],
            IupacBase::K => [0.0, 0.5, 0.5, 0.0],
            IupacBase::M => [0.5, 0.5, 0.0, 0.0],
            IupacBase::B => [0.0, 0.333, 0.333, 0.333],
            IupacBase::D => [0.333, 0.0, 0.333, 0.333],
            IupacBase::H => [0.333, 0.333, 0.0, 0.333],
            IupacBase::V => [0.333, 0.333, 0.333, 0.0],
            IupacBase::N => [0.25, 0.25, 0.25, 0.25],
        }
    }

    pub fn complement(&self) -> IupacBase {
        match self {
            IupacBase::A => IupacBase::T,
            IupacBase::C => IupacBase::G,
            IupacBase::G => IupacBase::C,
            IupacBase::T => IupacBase::A,
            IupacBase::R => IupacBase::Y,
            IupacBase::Y => IupacBase::R,
            IupacBase::S => IupacBase::S,
            IupacBase::W => IupacBase::W,
            IupacBase::K => IupacBase::M,
            IupacBase::M => IupacBase::K,
            IupacBase::B => IupacBase::V,
            IupacBase::D => IupacBase::H,
            IupacBase::H => IupacBase::D,
            IupacBase::V => IupacBase::B,
            IupacBase::N => IupacBase::N,
        }
    }

    pub fn to_string(&self) -> &'static str {
        match self {
            IupacBase::A => "A",
            IupacBase::C => "C",
            IupacBase::G => "G",
            IupacBase::T => "T",
            IupacBase::R => "R",
            IupacBase::Y => "Y",
            IupacBase::S => "S",
            IupacBase::W => "W",
            IupacBase::K => "K",
            IupacBase::M => "M",
            IupacBase::B => "B",
            IupacBase::D => "D",
            IupacBase::H => "H",
            IupacBase::V => "V",
            IupacBase::N => "N",
        }
    }

    pub fn from_char(c: char) -> Result<IupacBase, anyhow::Error> {
        match c {
            'A' => Ok(IupacBase::A),
            'C' => Ok(IupacBase::C),
            'G' => Ok(IupacBase::G),
            'T' => Ok(IupacBase::T),
            'R' => Ok(IupacBase::R),
            'Y' => Ok(IupacBase::Y),
            'S' => Ok(IupacBase::S),
            'W' => Ok(IupacBase::W),
            'K' => Ok(IupacBase::K),
            'M' => Ok(IupacBase::M),
            'B' => Ok(IupacBase::B),
            'D' => Ok(IupacBase::D),
            'H' => Ok(IupacBase::H),
            'V' => Ok(IupacBase::V),
            'N' => Ok(IupacBase::N),
            _ => bail!("Invalid IUPAC base: {}", c),
        }
    }

    pub fn contain_other_iupac_base(&self, other: IupacBase) -> bool {
        match self {
            IupacBase::A => match other {
                IupacBase::A => true,
                _ => false,
            },
            IupacBase::C => match other {
                IupacBase::C => true,
                _ => false,
            },
            IupacBase::G => match other {
                IupacBase::G => true,
                _ => false,
            },
            IupacBase::T => match other {
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::R => match other {
                IupacBase::A => true,
                IupacBase::G => true,
                _ => false,
            },
            IupacBase::Y => match other {
                IupacBase::C => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::S => match other {
                IupacBase::G => true,
                IupacBase::C => true,
                _ => false,
            },
            IupacBase::W => match other {
                IupacBase::A => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::K => match other {
                IupacBase::G => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::M => match other {
                IupacBase::A => true,
                IupacBase::C => true,
                _ => false,
            },
            IupacBase::B => match other {
                IupacBase::C => true,
                IupacBase::G => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::D => match other {
                IupacBase::A => true,
                IupacBase::G => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::H => match other {
                IupacBase::A => true,
                IupacBase::C => true,
                IupacBase::T => true,
                _ => false,
            },
            IupacBase::V => match other {
                IupacBase::A => true,
                IupacBase::C => true,
                IupacBase::G => true,
                _ => false,
            },
            IupacBase::N => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_regex() {
        for base in vec![
            IupacBase::A,
            IupacBase::C,
            IupacBase::G,
            IupacBase::T,
            IupacBase::R,
            IupacBase::Y,
            IupacBase::S,
            IupacBase::W,
            IupacBase::K,
            IupacBase::M,
            IupacBase::B,
            IupacBase::D,
            IupacBase::H,
            IupacBase::V,
            IupacBase::N,
        ] {
            let regex = base.to_regex();
            let expected = match base {
                IupacBase::A => "A",
                IupacBase::C => "C",
                IupacBase::G => "G",
                IupacBase::T => "T",
                IupacBase::R => "[AG]",
                IupacBase::Y => "[CT]",
                IupacBase::S => "[GC]",
                IupacBase::W => "[AT]",
                IupacBase::K => "[GT]",
                IupacBase::M => "[AC]",
                IupacBase::B => "[CGT]",
                IupacBase::D => "[AGT]",
                IupacBase::H => "[ACT]",
                IupacBase::V => "[ACG]",
                IupacBase::N => ".",
            };
            assert_eq!(regex, expected)
        }
    }

    #[test]
    fn test_complement() {
        for base in vec![
            IupacBase::A,
            IupacBase::C,
            IupacBase::G,
            IupacBase::T,
            IupacBase::R,
            IupacBase::Y,
            IupacBase::S,
            IupacBase::W,
            IupacBase::K,
            IupacBase::M,
            IupacBase::B,
            IupacBase::D,
            IupacBase::H,
            IupacBase::V,
            IupacBase::N,
        ] {
            let complement = base.complement();
            let expected = match base {
                IupacBase::A => IupacBase::T,
                IupacBase::C => IupacBase::G,
                IupacBase::G => IupacBase::C,
                IupacBase::T => IupacBase::A,
                IupacBase::R => IupacBase::Y,
                IupacBase::Y => IupacBase::R,
                IupacBase::S => IupacBase::S,
                IupacBase::W => IupacBase::W,
                IupacBase::K => IupacBase::M,
                IupacBase::M => IupacBase::K,
                IupacBase::B => IupacBase::V,
                IupacBase::D => IupacBase::H,
                IupacBase::H => IupacBase::D,
                IupacBase::V => IupacBase::B,
                IupacBase::N => IupacBase::N,
            };
            assert_eq!(complement, expected)
        }
    }
}
