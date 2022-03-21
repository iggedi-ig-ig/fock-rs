use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::window::Window;
use log::{info, LevelFilter};
use nalgebra::{Point3, Translation3, Vector3};
use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;
use scf::SelfConsistentField;
use std::sync::{Arc, Mutex};
use std::thread;

#[derive(Copy, Clone)]
struct DataPoint {
    position: Point3<f32>,
    color: Point3<f32>,
    prob: f32,
}

const POINTS_PER_N: usize = 2_500_000;
const POINTS_PER_ITER: usize = 125_000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::builder()
        .filter_level(LevelFilter::Debug)
        .init();

    let molecule = chemfiles::xyz::read_xyz_file(
        "../molecules/naphthalene.xyz",
        &basis_set::basis_sets::BASIS_STO_3G,
    )?;
    if let Some(result) = molecule.try_scf(100, 1e-6, 0) {
        let n_basis = result.n_basis;

        let mut rng = XorShiftRng::from_entropy();
        let mut window = Window::new("window");

        molecule.iter().for_each(|atom| {
            let mut sphere = window.add_sphere(0.25);
            let color = atom.atom_type().color();
            sphere.set_color(color[0], color[1], color[2]);
            sphere.set_local_translation(Translation3::from(atom.position().map(|f| f as f32)))
        });

        let data_points = Arc::new(Mutex::new(
            (0..n_basis).map(|_| Vec::new()).collect::<Vec<_>>(),
        ));
        let n = Arc::new(Mutex::new(0));

        let mut min_prob = 0.02;
        let (min, max): (Vector3<f64>, Vector3<f64>) =
            molecule
                .iter()
                .fold((Vector3::zeros(), Vector3::zeros()), |(min, max), curr| {
                    const EPSILON: f64 = 1e-4;
                    let max_r = curr
                        .basis()
                        .iter()
                        .flat_map(|basis| basis.primitives())
                        .fold(0.0f64, |f, curr| {
                            f.max(f64::sqrt(
                                (f64::ln(curr.coefficient()) - f64::ln(EPSILON)) / curr.exponent(),
                            ))
                        });
                    (
                        min.zip_map(&curr.position(), |a, b| a.min(b - max_r)),
                        max.zip_map(&curr.position(), |a, b| a.max(b + max_r)),
                    )
                });

        thread::spawn({
            let data_points = data_points.clone();
            let n = n.clone();

            move || {
                while let Some(n) = {
                    let lock = data_points.lock().unwrap();

                    let curr_n = n.lock().unwrap();
                    if lock[(*curr_n)].len() < POINTS_PER_N {
                        Some(*curr_n)
                    } else {
                        (0..n_basis).find(|n| lock[*n].len() < POINTS_PER_N)
                    }
                } {
                    for _ in 0..POINTS_PER_ITER {
                        let point = min + rng.gen::<Vector3<f64>>().component_mul(&(max - min));
                        let wave = result.orbitals[n].evaluate(&point);
                        let prob = wave.powi(2);

                        let mut lock = data_points.lock().unwrap();
                        lock[n].push(DataPoint {
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

                info!("All energy levels have been filled!");
            }
        });

        while window.render() {
            let current_energy_level = *n.lock().unwrap();
            window.set_title(&*format!(
                "Energy Level: {}/{} (E: {:+0.5} Hartrees) | total energy: {:+0.5} Hartrees",
                current_energy_level,
                n_basis - 1,
                result.orbital_energies[current_energy_level], // <-- Hartree to eV conversion factor,
                result.total_energy
            ));

            window.events().iter().for_each(|event| {
                if let WindowEvent::Key(key, Action::Press, _) = event.value {
                    match key {
                        Key::Up => min_prob *= 1.25,
                        Key::Down => min_prob /= 1.25,
                        Key::Left => {
                            *n.lock().unwrap() = {
                                if current_energy_level > 0 {
                                    current_energy_level - 1
                                } else {
                                    n_basis - 1
                                }
                            }
                        }
                        Key::Right => {
                            *n.lock().unwrap() = {
                                if current_energy_level < n_basis - 1 {
                                    current_energy_level + 1
                                } else {
                                    0
                                }
                            };
                        }
                        _ => {}
                    }
                }
            });

            let lock = data_points.lock().unwrap();
            lock[current_energy_level]
                .iter()
                .filter(|point| point.prob > min_prob)
                .for_each(|point| window.draw_point(&point.position, &point.color));
        }
        Ok(())
    } else {
        panic!("SCF didn't converge!");
    }
}
