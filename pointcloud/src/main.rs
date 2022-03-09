use crate::molecules::*;
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

const POINTS_PER_N: usize = 2_500_000;
const POINTS_PER_ITER: usize = 50_000;

fn main() {
    let basis_set = &basis_set::basis_sets::BASIS_6_31G;
    let molecule = WATER.build(basis_set);
    if let Some(result) = molecule.try_scf(5000, 1e-6, 0) {
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
                        Vector3::new(
                            min.x.min(curr.position().x - max_r),
                            min.y.min(curr.position().y - max_r),
                            min.z.min(curr.position().z - max_r),
                        ),
                        Vector3::new(
                            max.x.max(curr.position().x + max_r),
                            max.y.max(curr.position().y + max_r),
                            max.z.max(curr.position().z + max_r),
                        ),
                    )
                });
        while window.render() {
            window.set_title(&*format!(
                "Energy Level: {}/{} (E: {:+0.2}eV) | total energy: {:+0.2} Hartrees",
                n,
                n_basis - 1,
                result.orbital_energies[n] * 27.211, // <-- Hartree to eV conversion factor,
                result.total_energy
            ));
            let curr_n = if data_points[n].len() < POINTS_PER_N {
                Some(n)
            } else {
                (0..n_basis).find(|n| data_points[*n].len() < POINTS_PER_N)
            };
            if let Some(n) = curr_n {
                if data_points[n].len() < POINTS_PER_N {
                    for _ in 0..POINTS_PER_ITER {
                        let point = min + rng.gen::<Vector3<f64>>().component_mul(&(max - min));
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

                    window.draw_text(
                        &*format!(
                            "Currently filling energy level {}\n {} points",
                            n,
                            data_points[n].len()
                        ),
                        &Point2::new(10.0, 5.0),
                        window.width() as f32 / 25.0,
                        &Font::default(),
                        &Point3::new(1.0, 1.0, 1.0),
                    );
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
        }
    } else {
        panic!("SCF didn't converge!");
    }
}
