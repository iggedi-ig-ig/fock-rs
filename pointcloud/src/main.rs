use basis_set::periodic_table::AtomType::{Hydrogen, Oxygen};
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
    let basis = &basis_set::basis_sets::BASIS_STO_3G;
    let a = 104.5f64.to_radians() * 0.5;
    let l = 0.96 * 1.89;
    let molecule = [
        basis.get(Vector3::new(-a.sin(), -a.cos(), 0.0) * l, Hydrogen, 0),
        basis.get(Vector3::new(a.sin(), -a.cos(), 0.0) * l, Hydrogen, 0),
        basis.get(Vector3::zeros(), Oxygen, 0),
    ];
    if let Some(result) = molecule.try_scf(1000, 1e-12) {
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
            let [r, g, b] = atom.atom_type().color();
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
                    "point cloud size: {}\nenergy level: {}/{} (E: {:0.4} eV)\nmin prob: {:0.5}",
                    points.len(),
                    n,
                    result.wave_function.basis_functions().len() - 1,
                    result.orbital_energies[n] * 27.211,
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
