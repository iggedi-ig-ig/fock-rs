use serde::Deserialize;
use std::mem::transmute;

#[repr(usize)]
#[derive(Deserialize, Debug, Eq, PartialEq, Hash, Copy, Clone)]
pub enum AtomType {
    #[serde(rename = "1")]
    Hydrogen = 1,
    #[serde(rename = "2")]
    Helium = 2,
    #[serde(rename = "3")]
    Lithium = 3,
    #[serde(rename = "4")]
    Beryllium = 4,
    #[serde(rename = "5")]
    Boron = 5,
    #[serde(rename = "6")]
    Carbon = 6,
    #[serde(rename = "7")]
    Nitrogen = 7,
    #[serde(rename = "8")]
    Oxygen = 8,
    #[serde(rename = "9")]
    Fluorine = 9,
    #[serde(rename = "10")]
    Neon = 10,
    #[serde(rename = "11")]
    Sodium = 11,
    #[serde(rename = "12")]
    Magnesium = 12,
    #[serde(rename = "13")]
    Aluminium = 13,
    #[serde(rename = "14")]
    Silicon = 14,
    #[serde(rename = "15")]
    Phosphorous = 15,
    #[serde(rename = "16")]
    Sulfur = 16,
    #[serde(rename = "17")]
    Chlorine = 17,
    #[serde(rename = "18")]
    Argon = 18,
    #[serde(rename = "19")]
    Potassium = 19,
    #[serde(rename = "20")]
    Calcium = 20,
    #[serde(rename = "21")]
    Scandium = 21,
    #[serde(rename = "22")]
    Titanium = 22,
    #[serde(rename = "23")]
    Vanadium = 23,
    #[serde(rename = "24")]
    Chromium = 24,
    #[serde(rename = "25")]
    Manganese = 25,
    #[serde(rename = "26")]
    Iron = 26,
    #[serde(rename = "27")]
    Cobalt = 27,
    #[serde(rename = "28")]
    Nickel = 28,
    #[serde(rename = "29")]
    Copper = 29,
    #[serde(rename = "30")]
    Zinc = 30,
    #[serde(rename = "31")]
    Gallium = 31,
    #[serde(rename = "32")]
    Germanium = 32,
    #[serde(rename = "33")]
    Arsenic = 33,
    #[serde(rename = "34")]
    Selenium = 34,
    #[serde(rename = "35")]
    Bromine = 35,
    #[serde(rename = "36")]
    Krypton = 36,
}

impl AtomType {
    pub fn from_ordinal(ordinal: usize) -> Self {
        unsafe { transmute(ordinal) }
    }
}

#[test]
pub fn test_from_ordinal() {
    use crate::periodic_table::AtomType::*;

    assert_eq!(Hydrogen, AtomType::from_ordinal(1));
    assert_eq!(Oxygen, AtomType::from_ordinal(8));
}
