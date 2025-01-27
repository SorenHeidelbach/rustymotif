use anyhow::bail;
use std::fmt::Display;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
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
