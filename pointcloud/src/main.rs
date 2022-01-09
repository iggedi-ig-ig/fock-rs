use crate::molecules::{AMMONIA, BENZENE, METHANE, NITRITE, WATER};
use basis_set::periodic_table::AtomType;
use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::text::Font;
use kiss3d::window::Window;
use nalgebra::{Point2, Point3, Translation3, Vector3};
use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;
use scf::SelfConsistentField;

mod molecules;

#[derive(Copy, Clone)]
struct DataPoint {
    position: Point3<f32>,
    color: Point3<f32>,
    prob: f32,
}

fn main() {
    let basis_set = &basis_set::basis_sets::BASIS_3_21G;
    let molecule = [
        basis_set.get(Vector3::zeros(), AtomType::Hydrogen),
        basis_set.get(Vector3::new(bohr!(1.25), 0.0, 0.0), AtomType::Chlorine),
    ];
    if let Some(result) = molecule.try_scf(1000, 1e-6, 0) {
        let n_basis = result.orbitals.basis_functions().len();

        let mut rng = XorShiftRng::from_entropy();
        let mut window = Window::new("window");

        molecule.iter().for_each(|atom| {
            let mut sphere = window.add_sphere(0.25);
            let color = atom.atom_type().color();
            sphere.set_color(color[0], color[1], color[2]);
            sphere.set_local_translation(Translation3::from(atom.position().map(|f| f as f32)))
        });

        let mut data_points = (0..n_basis).map(|_| Vec::new()).collect::<Vec<_>>();
        let mut n = 0;
        let mut min_prob = 0.02;
        while window.render() {
            window.set_title(&*format!(
                "Energy Level: {}/{} (E: {:+0.3}eV) | total energy: {:+0.3} Hartrees",
                n,
                n_basis - 1,
                result.orbital_energies[n] * 27.211, // <-- Hartree to eV conversion factor,
                result.total_energy
            ));
            if data_points[n].len() < 1_000_000 {
                for _ in 0..50_000 {
                    let point = (rng.gen::<Vector3<f64>>() - Vector3::repeat(0.5))
                        .component_mul(&Vector3::new(7.5, 5.0, 7.5))
                        * 2.0;
                    let wave = result.orbitals.evaluate(point, n);
                    let prob = wave.powi(2);

                    data_points[n].push(DataPoint {
                        position: Point3::from(point.map(|f| f as f32)),
                        color: if wave > 0.0 {
                            Point3::new(1.0, 0.2, 0.2)
                        } else {
                            Point3::new(0.2, 0.2, 1.0)
                        },
                        prob: prob as f32,
                    });
                }
            }

            window.events().iter().for_each(|event| {
                if let WindowEvent::Key(key, Action::Press, _) = event.value {
                    match key {
                        Key::Up => min_prob *= 1.25,
                        Key::Down => min_prob /= 1.25,
                        Key::Left => n = if n > 0 { n - 1 } else { n_basis - 1 },
                        Key::Right => n = if n < n_basis - 1 { n + 1 } else { 0 },
                        _ => {}
                    }
                }
            });

            data_points[n]
                .iter()
                .filter(|point| point.prob > min_prob)
                .for_each(|point| window.draw_point(&point.position, &point.color));

            window.draw_text(
                &*format!("{} points", data_points[n].len()),
                &Point2::new(10.0, 5.0),
                window.width() as f32 / 25.0,
                &Font::default(),
                &Point3::new(1.0, 1.0, 1.0),
            );
        }
    } else {
        panic!("SCF didn't converge!");
    }
}
