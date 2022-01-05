use crate::BasisSet;

pub const DATA_STO3G: &str = include_str!("data/STO-3G.json");
pub const DATA_3_21G: &str = include_str!("data/3-21G.json");
pub const DATA_6_31G: &str = include_str!("data/6-31G.json");
pub const DATA_6_31G_ST: &str = include_str!("data/6-31G_st.json");

lazy_static::lazy_static! {
    pub static ref BASIS_STO_3G: BasisSet = serde_json::from_str(DATA_STO3G).expect("failed to read STO-3G basis set");
    pub static ref BASIS_3_21G: BasisSet = serde_json::from_str(DATA_3_21G).expect("failed to read 3-21G");
    pub static ref BASIS_6_31G: BasisSet = serde_json::from_str(DATA_6_31G).expect("failed to read 6-31G basis set");
    pub static ref BASIS_6_31G_ST: BasisSet = serde_json::from_str(DATA_6_31G_ST).expect("failed to read 6-31G* basis set");
}
