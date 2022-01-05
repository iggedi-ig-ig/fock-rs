use basis_set::basis_sets::{BASIS_6_31G, BASIS_STO_3G};
use basis_set::periodic_table::AtomType;
use kiss3d::camera::ArcBall;
use kiss3d::event::{Action, Key, WindowEvent};
use kiss3d::light::Light;
use kiss3d::text::Font;
use kiss3d::window::Window;
use nalgebra::{Point2, Point3, Translation3, Vector3};
use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;
use scf::SelfConsistentField;

const POINT_CLOUD_SIZE: usize = 50_000;
const POINT_CLOUD_ITER_SIZE: usize = 75_000;
const RENDER_SCALE: f32 = 1.0;

fn main() {
    let basis = &basis_set::basis_sets::BASIS_6_31G;
    let molecule = [
        basis
            .get(Vector3::new(-0.7, 0.0, 0.0), AtomType::Helium, 0)
            .unwrap(),
        basis
            .get(Vector3::new(0.7, 0.0, 0.0), AtomType::Hydrogen, 0)
            .unwrap(),
    ];
    if let Some(result) = molecule.try_scf(2, 1000, 1e-12) {
        println!("SCF-Cycle converged after {} iterations", result.iterations);
        print!(
            "Orbital energies: {:0.4}",
            result.orbital_energies.transpose()
        );
        print!("Density matrix: {:0.4}", result.density_matrix);
        print!(
            "Orbital matrix: {:0.4}",
            result.wave_function.coeff_matrix()
        );

        let mut window = Window::new("test");
        let mut rng = XorShiftRng::from_entropy();
        window.set_light(Light::StickToCamera);

        let fov = std::f32::consts::PI / 4.0;
        let mut cam = ArcBall::new_with_frustrum(
            fov,
            0.1,
            1024.0,
            Point3::new(0.0f32, 0.0, -1.0),
            Point3::origin(),
        );

        molecule.iter().for_each(|atom| {
            let mut sphere = window.add_sphere(0.25);
            // TODO
            let [r, g, b] = [1.0, 1.0, 1.0];
            sphere.set_color(r, g, b);
            sphere.set_local_translation(Translation3::from(
                Vector3::new(atom.position().x, atom.position().y, atom.position().z)
                    .map(|f| f as f32 * RENDER_SCALE),
            ));
        });
        let n_basis = result.wave_function.basis_functions().len();

        let mut n = 0;
        let mut points = Vec::new();
        let mut min_prob = 0.01;

        while window.render_with_camera(&mut cam) {
            if points.len() < POINT_CLOUD_SIZE {
                (0..POINT_CLOUD_ITER_SIZE).for_each(|_| {
                    let point = (rng.gen::<Vector3<f64>>() - Vector3::repeat(0.5)) * 2.0 * 5.0;
                    let wave = result.wave_function.evaluate(point, n);
                    let prob = wave.powi(2);

                    if prob > min_prob {
                        points.push((
                            Point3::from(point.map(|f| f as f32 * RENDER_SCALE)),
                            if wave > 0.0 {
                                Point3::new(0.2, 0.2, 1.0)
                            } else {
                                Point3::new(1.0, 0.2, 0.2)
                            },
                        ))
                    }
                });
            }

            points
                .iter()
                .for_each(|(point, color)| window.draw_point(point, color));

            window.draw_text(
                &*format!(
                    "point cloud size: {}\nenergy level: {}/{} (E: {:0.4} Hartrees)\nmin prob: {:0.5}",
                    points.len(),
                    n,
                    result.wave_function.basis_functions().len() - 1,
                    result.orbital_energies[n],
                    min_prob,
                ),
                &Point2::new(10.0, 10.0),
                50.0,
                &Font::default(),
                &Point3::new(1.0, 1.0, 1.0),
            );

            for x in window.events().iter() {
                if let WindowEvent::Key(key, Action::Press, _) = x.value {
                    match key {
                        Key::Left => {
                            n = if n > 0 { n - 1 } else { n_basis - 1 };
                            points.clear();
                        }
                        Key::Right => {
                            n = if n < n_basis - 1 { n + 1 } else { 0 };
                            points.clear();
                        }
                        Key::Up => {
                            min_prob *= 1.25;
                            points.clear();
                        }
                        Key::Down => {
                            min_prob *= 0.8;
                            points.clear();
                        }
                        Key::R => {}
                        _ => {}
                    }
                }
            }
        }
    }
}
