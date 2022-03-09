macro_rules! define_basis_sets {
    ($($name:ident => $path:literal),*) => {
        lazy_static::lazy_static! {
            $(
                pub static ref $name: crate::BasisSet =
                    serde_json::from_str(include_str!($path)).expect(&format!("Failed to read {}", $path));
            )*
        }
    };
}

define_basis_sets! {
    BASIS_STO_3G => "data/STO-3G.json",
    BASIS_3_21G => "data/3-21G.json",
    BASIS_6_31G => "data/6-31G.json",
    BASIS_6_311G => "data/6-311G.json",
    BASIS_6_31_PP_G => "data/6-31++G.json",
    BASIS_6_31G_ST => "data/6-31G_st.json"
}
