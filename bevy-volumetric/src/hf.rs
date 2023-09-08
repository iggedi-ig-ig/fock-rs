use bevy::{
    prelude::*,
    tasks::{AsyncComputeTaskPool, Task},
};
use futures_lite::future;
use scf::{HartreeFockResult, SelfConsistentField};

use crate::molecule::LoadedMolecule;

pub struct HfPlugin;

#[derive(Resource, Default, Clone)]
pub enum ConvergedScf {
    #[default]
    Unsolved,
    Solving,
    Solved(HartreeFockResult),
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, States)]
pub enum ScfConvergedState {
    Solving,
    Converged,
    #[default]
    Failed,
}

#[derive(Resource, Copy, Clone, Debug)]
pub struct ScfSettings {
    max_iters: usize,
    epsilon: f64,
    molecule_charge: i32,
}

impl ScfSettings {
    pub fn ion(molecule_charge: i32) -> Self {
        Self {
            molecule_charge,
            ..Default::default()
        }
    }
}

#[derive(Component, Deref, DerefMut)]
struct ScfTask(Task<Option<HartreeFockResult>>);

impl Default for ScfSettings {
    fn default() -> Self {
        Self {
            max_iters: 100,
            epsilon: 1e-6,
            molecule_charge: 0,
        }
    }
}

impl Plugin for HfPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Update, (solve_scf, poll_scf_tasks))
            .init_resource::<ScfSettings>()
            .init_resource::<ConvergedScf>()
            .add_state::<ScfConvergedState>();
    }
}

fn solve_scf(
    mut commands: Commands,
    mut solved_state: ResMut<NextState<ScfConvergedState>>,
    mut converged: ResMut<ConvergedScf>,
    molecule: Res<LoadedMolecule>,
    settings: Res<ScfSettings>,
) {
    let pool = AsyncComputeTaskPool::get();

    if molecule.is_changed() {
        let Some(molecule) = &**molecule else {return;};
        solved_state.set(ScfConvergedState::Solving);
        *converged = ConvergedScf::Solving;

        let molecule = molecule.clone();
        let settings = *settings;
        commands.spawn(ScfTask(pool.spawn(async move {
            molecule.try_scf(
                settings.max_iters,
                settings.epsilon,
                settings.molecule_charge,
            )
        })));
    }
}

fn poll_scf_tasks(
    mut commands: Commands,
    mut tasks: Query<(Entity, &mut ScfTask)>,
    mut solved_state: ResMut<NextState<ScfConvergedState>>,
    mut converged: ResMut<ConvergedScf>,
) {
    for (entity, mut task) in &mut tasks {
        if let Some(result) = future::block_on(future::poll_once(&mut task.0)) {
            info!("done computing scf");
            match result {
                Some(result) => {
                    solved_state.set(ScfConvergedState::Converged);
                    *converged = ConvergedScf::Solved(result);
                }
                None => {
                    solved_state.set(ScfConvergedState::Failed);
                    *converged = ConvergedScf::Unsolved;
                }
            }

            commands.entity(entity).despawn_recursive();
        }
    }
}
