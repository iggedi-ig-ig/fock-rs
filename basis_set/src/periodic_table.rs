use crate::AtomType::*;
use serde::Deserialize;
use std::mem::transmute;

// TODO: there has to be a cleaner way to do this
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
    #[serde(rename = "37")]
    Rubidium = 37,
    #[serde(rename = "38")]
    Strontium = 38,
    #[serde(rename = "39")]
    Yttrium = 39,
    #[serde(rename = "40")]
    Zirconium = 40,
    #[serde(rename = "41")]
    Niobium = 41,
    #[serde(rename = "42")]
    Molybdenum = 42,
    #[serde(rename = "43")]
    Technetium = 43,
    #[serde(rename = "44")]
    Ruthenium = 44,
    #[serde(rename = "45")]
    Rhodium = 45,
    #[serde(rename = "46")]
    Palladium = 46,
    #[serde(rename = "47")]
    Silver = 47,
    #[serde(rename = "48")]
    Cadmium = 48,
    #[serde(rename = "49")]
    Indium = 49,
    #[serde(rename = "50")]
    Tin = 50,
    #[serde(rename = "51")]
    Antimony = 51,
    #[serde(rename = "52")]
    Tellurium = 52,
    #[serde(rename = "53")]
    Iodine = 53,
    #[serde(rename = "54")]
    Xenon = 54,
    #[serde(rename = "55")]
    Caesium = 55,
    #[serde(rename = "56")]
    Barium = 56,
    #[serde(rename = "57")]
    Lanthanum = 57,
    #[serde(rename = "58")]
    Cerium = 58,
    #[serde(rename = "59")]
    Praseodymium = 59,
    #[serde(rename = "60")]
    Neodymium = 60,
    #[serde(rename = "61")]
    Promethium = 61,
    #[serde(rename = "62")]
    Samarium = 62,
    #[serde(rename = "63")]
    Europium = 63,
    #[serde(rename = "64")]
    Gadolinium = 64,
    #[serde(rename = "65")]
    Terbium = 65,
    #[serde(rename = "66")]
    Dysprosium = 66,
    #[serde(rename = "67")]
    Holmium = 67,
    #[serde(rename = "68")]
    Erbium = 68,
    #[serde(rename = "69")]
    Thulium = 69,
    #[serde(rename = "70")]
    Ytterbium = 70,
    #[serde(rename = "71")]
    Lutetium = 71,
    #[serde(rename = "72")]
    Hafnium = 72,
    #[serde(rename = "73")]
    Tantalum = 73,
    #[serde(rename = "74")]
    Tungsten = 74,
    #[serde(rename = "75")]
    Rhenium = 75,
    #[serde(rename = "76")]
    Osmium = 76,
    #[serde(rename = "77")]
    Iridium = 77,
    #[serde(rename = "78")]
    Platinum = 78,
    #[serde(rename = "79")]
    Gold = 79,
    #[serde(rename = "80")]
    Mercury = 80,
    #[serde(rename = "81")]
    Thallium = 81,
    #[serde(rename = "82")]
    Lead = 82,
    #[serde(rename = "83")]
    Bismuth = 83,
    #[serde(rename = "84")]
    Polonium = 84,
    #[serde(rename = "85")]
    Astatine = 85,
    #[serde(rename = "86")]
    Radon = 86,
    #[serde(rename = "87")]
    Francium = 87,
    #[serde(rename = "88")]
    Radium = 88,
    #[serde(rename = "89")]
    Actinium = 89,
    #[serde(rename = "90")]
    Thorium = 90,
    #[serde(rename = "91")]
    Protactinium = 91,
    #[serde(rename = "92")]
    Uranium = 92,
    #[serde(rename = "93")]
    Neptunium = 93,
    #[serde(rename = "94")]
    Plutonium = 94,
    #[serde(rename = "95")]
    Americium = 95,
    #[serde(rename = "96")]
    Curium = 96,
    #[serde(rename = "97")]
    Berkelium = 97,
    #[serde(rename = "98")]
    Californium = 98,
    #[serde(rename = "99")]
    Einsteinium = 99,
    #[serde(rename = "100")]
    Fermium = 100,
    #[serde(rename = "101")]
    Mendelevium = 101,
    #[serde(rename = "102")]
    Nobelium = 102,
    #[serde(rename = "103")]
    Lawrencium = 103,
    #[serde(rename = "104")]
    Rutherfordium = 104,
    #[serde(rename = "105")]
    Dubnium = 105,
    #[serde(rename = "106")]
    Seaborgium = 106,
    #[serde(rename = "107")]
    Bohrium = 107,
    #[serde(rename = "108")]
    Hassium = 108,
    #[serde(rename = "109")]
    Meitnerium = 109,
    #[serde(rename = "110")]
    Darmstadtium = 110,
    #[serde(rename = "111")]
    Roentgenium = 111,
    #[serde(rename = "112")]
    Copernicium = 112,
    #[serde(rename = "113")]
    Nohinium = 113,
    #[serde(rename = "114")]
    Flerovium = 114,
    #[serde(rename = "115")]
    Moscovium = 115,
    #[serde(rename = "116")]
    Livermorium = 116,
    #[serde(rename = "117")]
    Tennessium = 117,
    #[serde(rename = "118")]
    Oganesson = 118,
    #[serde(rename = "119")]
    Ununennium = 119,
    #[serde(rename = "120")]
    Unbinilium = 120,
    #[serde(rename = "121")]
    Unbiunium = 121,
    #[serde(rename = "122")]
    Unbibium = 122,
    #[serde(rename = "123")]
    Unbitrium = 123,
    #[serde(rename = "124")]
    Unbiquadium = 124,
    #[serde(rename = "125")]
    Unbipentium = 125,
    #[serde(rename = "126")]
    Unbihexium = 126,
    #[serde(rename = "127")]
    Unbiseptium = 127,
}

impl AtomType {
    pub fn from_ordinal(ordinal: usize) -> Self {
        unsafe { transmute(ordinal) }
    }

    pub fn color(&self) -> [f32; 3] {
        match *self {
            Hydrogen => [0.9; 3],
            Carbon => [0.1; 3],
            Oxygen => [0.9, 0.1, 0.1],
            Chlorine => [0.1, 0.9, 0.1],
            Nitrogen => [0.1, 0.1, 0.9],
            Fluorine => [188.0 / 255.0, 166.0 / 255.0, 65.0 / 255.0],
            Iron => [0.8; 3],
            _ => [0.5; 3],
        }
    }

    pub fn from_symbol(symbol: &str) -> Self {
        match symbol {
            "H" => Hydrogen,
            "He" => Helium,
            "C" => Carbon,
            "N" => Nitrogen,
            "O" => Oxygen,
            "F" => Fluorine,
            "Fe" => Iron,
            _ => panic!("not implemented yet"),
        }
    }
}

#[test]
pub fn test_from_ordinal() {
    use crate::periodic_table::AtomType::*;

    assert_eq!(Hydrogen, AtomType::from_ordinal(1));
    assert_eq!(Oxygen, AtomType::from_ordinal(8));
}
