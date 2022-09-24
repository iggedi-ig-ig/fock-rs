use bytemuck::{Pod, Zeroable};
use dolly::glam::{EulerRot, Vec3};
use dolly::prelude::*;
use log::{error, info, LevelFilter};
use nalgebra::{Vector2, Vector3};
use rayon::prelude::*;
use scf::molecular_orbitals::MolecularOrbitals;
use scf::{HartreeFockResult, SelfConsistentField};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::mem::size_of;
use std::num::NonZeroU32;
use std::time::Instant;
use wgpu::{
    include_spirv, AddressMode, Backends, BindGroup, BindGroupDescriptor, BindGroupEntry,
    BindGroupLayoutDescriptor, BindGroupLayoutEntry, BindingResource, BindingType, BlendState,
    Color, ColorTargetState, ColorWrites, CommandEncoderDescriptor, Device, DeviceDescriptor,
    Extent3d, Features, FilterMode, FragmentState, ImageDataLayout, Instance, Limits, LoadOp,
    Operations, PipelineLayoutDescriptor, PresentMode, PrimitiveState, PrimitiveTopology,
    PushConstantRange, Queue, RenderPassColorAttachment, RenderPassDescriptor, RenderPipeline,
    RenderPipelineDescriptor, RequestAdapterOptions, Sampler, SamplerBindingType,
    SamplerBorderColor, SamplerDescriptor, ShaderStages, Surface, SurfaceConfiguration,
    SurfaceError, Texture, TextureDescriptor, TextureDimension, TextureFormat, TextureSampleType,
    TextureUsages, TextureView, TextureViewDescriptor, TextureViewDimension, VertexState,
};
use winit::dpi::PhysicalSize;
use winit::event::{
    DeviceEvent, ElementState, Event, KeyboardInput, MouseScrollDelta, VirtualKeyCode, WindowEvent,
};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::window::{CursorGrabMode, Window, WindowBuilder};

#[repr(C)]
#[derive(Copy, Clone, Zeroable, Pod)]
struct PushConstants {
    dens_multiplier: f32,
    pos_x: f32,
    pos_y: f32,
    pos_z: f32,
    yaw: f32,
    pitch: f32,
    inv_aspect: f32,
    vox_side: u32,
}

const VOX_COUNT_SIDE: usize = 150;
const VOX_TEXTURE_SCALE: f64 = 20.0;

struct State {
    surface: Surface,
    device: Device,
    queue: Queue,
    config: SurfaceConfiguration,
    size: PhysicalSize<u32>,

    density_texture: Texture,

    bind_group: BindGroup,
    render_pipeline: RenderPipeline,

    last_time: Instant,

    mouse_movement: Vector2<f32>,
    keys: HashMap<VirtualKeyCode, bool>,

    cam_rig: CameraRig,
    density_scale: f32,

    hf_result: HartreeFockResult,
    energy_level: usize,
}

impl State {
    fn create_wave_texture_data(n: usize, orbitals: &MolecularOrbitals) -> Vec<f32> {
        (0..VOX_COUNT_SIDE.pow(3))
            .into_par_iter()
            .map(|i| {
                let x = i / VOX_COUNT_SIDE.pow(2);
                let y = (i / VOX_COUNT_SIDE) % VOX_COUNT_SIDE;
                let z = i % VOX_COUNT_SIDE;
                orbitals[n].evaluate(
                    &(Vector3::new(x, y, z).map(|x| x as f64 / VOX_COUNT_SIDE as f64 - 0.5)
                        * VOX_TEXTURE_SCALE),
                ) as f32
            })
            .collect()
    }

    fn update_density_texture(&mut self) {
        let before = Instant::now();

        self.queue.write_texture(
            self.density_texture.as_image_copy(),
            bytemuck::cast_slice(&Self::create_wave_texture_data(
                self.energy_level,
                &self.hf_result.orbitals,
            )),
            ImageDataLayout {
                offset: 0,
                bytes_per_row: NonZeroU32::new(VOX_COUNT_SIDE as u32 * size_of::<f32>() as u32),
                rows_per_image: NonZeroU32::new(VOX_COUNT_SIDE as u32),
            },
            Extent3d {
                width: VOX_COUNT_SIDE as u32,
                height: VOX_COUNT_SIDE as u32,
                depth_or_array_layers: VOX_COUNT_SIDE as u32,
            },
        );

        info!(
            "Energy level {}. Took {:?}",
            self.energy_level,
            before.elapsed()
        );
    }

    async fn new(window: &Window) -> Self {
        let size = window.inner_size();

        let molecule = chemfiles::xyz::read_xyz_file(
            "chemfiles/molecules/benzene.xyz",
            &basis_set::basis_sets::BASIS_STO_3G,
        )
        .expect("xyz file is invalid");

        let hf_result = molecule.try_scf(100, 1e-6, 0).expect("scf failed");

        let instance = Instance::new(Backends::VULKAN);
        let surface = unsafe { instance.create_surface(window) };
        let adapter = instance
            .request_adapter(&RequestAdapterOptions {
                power_preference: Default::default(),
                force_fallback_adapter: false,
                compatible_surface: Some(&surface),
            })
            .await
            .unwrap();

        info!("got adapter. GPU Name: {}", adapter.get_info().name);

        let (device, queue) = adapter
            .request_device(
                &DeviceDescriptor {
                    label: None,
                    features: Features::TEXTURE_ADAPTER_SPECIFIC_FORMAT_FEATURES
                        | Features::ADDRESS_MODE_CLAMP_TO_BORDER
                        | Features::PUSH_CONSTANTS,
                    limits: Limits {
                        max_push_constant_size: size_of::<PushConstants>() as u32,
                        ..Default::default()
                    },
                },
                None,
            )
            .await
            .unwrap();

        let supported_formats = surface.get_supported_formats(&adapter);

        info!(
            "{} supported formats: {:?}",
            supported_formats.len(),
            supported_formats
        );

        let config = SurfaceConfiguration {
            usage: TextureUsages::RENDER_ATTACHMENT,
            format: supported_formats[0],
            width: size.width,
            height: size.height,
            present_mode: PresentMode::Fifo,
        };
        surface.configure(&device, &config);

        let density_texture = device.create_texture(&TextureDescriptor {
            label: None,
            size: Extent3d {
                width: VOX_COUNT_SIDE as u32,
                height: VOX_COUNT_SIDE as u32,
                depth_or_array_layers: VOX_COUNT_SIDE as u32,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: TextureDimension::D3,
            format: TextureFormat::R32Float,
            usage: TextureUsages::COPY_DST | TextureUsages::TEXTURE_BINDING,
        });
        let density_view = density_texture.create_view(&TextureViewDescriptor::default());
        let density_sampler = device.create_sampler(&SamplerDescriptor {
            address_mode_u: AddressMode::ClampToBorder,
            address_mode_v: AddressMode::ClampToBorder,
            address_mode_w: AddressMode::ClampToBorder,
            mag_filter: FilterMode::Linear,
            min_filter: FilterMode::Linear,
            mipmap_filter: FilterMode::Linear,
            border_color: Some(SamplerBorderColor::TransparentBlack),
            ..Default::default()
        });

        let bg_layout = device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            label: None,
            entries: &[
                BindGroupLayoutEntry {
                    binding: 0,
                    visibility: ShaderStages::FRAGMENT,
                    ty: BindingType::Texture {
                        sample_type: TextureSampleType::Float { filterable: true },
                        view_dimension: TextureViewDimension::D3,
                        multisampled: false,
                    },
                    count: None,
                },
                BindGroupLayoutEntry {
                    binding: 1,
                    visibility: ShaderStages::FRAGMENT,
                    ty: BindingType::Sampler(SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });
        let bind_group = device.create_bind_group(&BindGroupDescriptor {
            label: None,
            layout: &bg_layout,
            entries: &[
                BindGroupEntry {
                    binding: 0,
                    resource: BindingResource::TextureView(&density_view),
                },
                BindGroupEntry {
                    binding: 1,
                    resource: BindingResource::Sampler(&density_sampler),
                },
            ],
        });

        let vert_shader = device.create_shader_module(include_spirv!("../shaders/vert.spv"));
        let frag_shader = device.create_shader_module(include_spirv!("../shaders/frag.spv"));

        let rp_layout = device.create_pipeline_layout(&PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[&bg_layout],
            push_constant_ranges: &[PushConstantRange {
                stages: ShaderStages::FRAGMENT,
                range: 0..size_of::<PushConstants>() as u32,
            }],
        });
        let render_pipeline = device.create_render_pipeline(&RenderPipelineDescriptor {
            label: None,
            layout: Some(&rp_layout),
            vertex: VertexState {
                module: &vert_shader,
                entry_point: "main",
                buffers: &[],
            },
            primitive: PrimitiveState {
                topology: PrimitiveTopology::TriangleStrip,
                ..Default::default()
            },
            depth_stencil: None,
            multisample: Default::default(),
            fragment: Some(FragmentState {
                module: &frag_shader,
                entry_point: "main",
                targets: &[Some(ColorTargetState {
                    format: config.format,
                    blend: Some(BlendState::REPLACE),
                    write_mask: ColorWrites::ALL,
                })],
            }),
            multiview: None,
        });

        Self {
            surface,
            device,
            queue,
            config,
            size,

            bind_group,
            render_pipeline,

            density_texture,

            last_time: Instant::now(),

            mouse_movement: Vector2::zeros(),
            keys: HashMap::default(),

            cam_rig: CameraRig::builder()
                .with(Position::new(Vec3::Y))
                .with(YawPitch::new())
                .with(Smooth::new_position(0.1))
                .build(),
            density_scale: 1.0,

            hf_result,
            energy_level: 0,
        }
    }

    fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.size = new_size;
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
        }
    }

    fn input(&mut self, event: &WindowEvent) {
        if let WindowEvent::KeyboardInput {
            input:
                KeyboardInput {
                    state,
                    virtual_keycode: Some(key),
                    ..
                },
            ..
        } = event
        {
            match key {
                VirtualKeyCode::Left => {
                    if let ElementState::Pressed = state {
                        if self.energy_level > 0 {
                            self.energy_level -= 1;
                        } else {
                            self.energy_level = self.hf_result.n_basis - 1;
                        }

                        self.update_density_texture();
                    }
                }
                VirtualKeyCode::Right => {
                    if let ElementState::Pressed = state {
                        if self.energy_level < self.hf_result.n_basis - 1 {
                            self.energy_level += 1;
                        } else {
                            self.energy_level = 0;
                        }

                        self.update_density_texture();
                    }
                }
                key => {
                    let pressed = match state {
                        ElementState::Pressed => true,
                        ElementState::Released => false,
                    };

                    self.keys
                        .entry(*key)
                        .and_modify(|entry| *entry = pressed)
                        .or_insert(pressed);
                    return;
                }
            }
        } else if let &WindowEvent::MouseWheel {
            delta: MouseScrollDelta::LineDelta(_x, y),
            ..
        } = event
        {
            match y.total_cmp(&0.0) {
                Ordering::Less => self.density_scale /= 1.25,
                Ordering::Equal => {}
                Ordering::Greater => self.density_scale *= 1.25,
            }
        }
    }

    fn update(&mut self) {
        // update routine

        let is_down = |key| *self.keys.get(&key).unwrap_or(&false);
        let move_axis =
            |k1, k2| if is_down(k1) { 1.0 } else { 0.0 } - if is_down(k2) { 1.0 } else { 0.0 };

        let move_fwd = move_axis(VirtualKeyCode::W, VirtualKeyCode::S);
        let move_right = move_axis(VirtualKeyCode::D, VirtualKeyCode::A);
        let move_up = move_axis(VirtualKeyCode::E, VirtualKeyCode::Q);

        let move_vec = self.cam_rig.final_transform.rotation
            * Vec3::new(move_right, move_up, move_fwd).clamp_length_max(1.0)
            * if is_down(VirtualKeyCode::LShift) {
                0.25
            } else {
                0.05
            };

        self.cam_rig
            .driver_mut::<YawPitch>()
            .rotate_yaw_pitch(self.mouse_movement.x, self.mouse_movement.y);
        self.cam_rig.driver_mut::<Position>().translate(move_vec);
        self.cam_rig.update(self.last_time.elapsed().as_secs_f32());
        self.last_time = Instant::now();

        self.mouse_movement *= 0.0;
    }

    fn render(&mut self) -> Result<(), SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&CommandEncoderDescriptor {
                label: Some("Command Encoder"),
            });
        {
            let pos_dolly = self.cam_rig.final_transform.position;
            let (yaw, pitch, _) = self
                .cam_rig
                .final_transform
                .rotation
                .to_euler(EulerRot::YXZ);

            let mut render_pass = encoder.begin_render_pass(&RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: Operations {
                        load: LoadOp::Clear(Color::BLACK),
                        store: true,
                    },
                })],
                depth_stencil_attachment: None,
            });

            render_pass.set_pipeline(&self.render_pipeline);
            render_pass.set_push_constants(
                ShaderStages::FRAGMENT,
                0,
                bytemuck::bytes_of(&PushConstants {
                    dens_multiplier: self.density_scale,
                    pos_x: pos_dolly.x,
                    pos_y: pos_dolly.y,
                    pos_z: pos_dolly.z,
                    yaw,
                    pitch,
                    inv_aspect: self.size.height as f32 / self.size.width as f32,
                    vox_side: VOX_COUNT_SIDE as u32,
                }),
            );
            render_pass.set_bind_group(0, &self.bind_group, &[]);
            render_pass.draw(0..4, 0..1);
        }

        self.queue.submit(Some(encoder.finish()));
        output.present();
        Ok(())
    }
}

#[tokio::main]
async fn main() {
    pretty_env_logger::formatted_builder()
        .filter_level(LevelFilter::Info)
        .init();

    let event_loop = EventLoop::new();
    let window = WindowBuilder::new()
        .with_inner_size(PhysicalSize::new(1920, 1080))
        .build(&event_loop)
        .unwrap();

    let mut state = State::new(&window).await;
    state.update_density_texture();

    window
        .set_cursor_grab(CursorGrabMode::Confined)
        .expect("failed to grab cursor");
    window.set_cursor_visible(false);

    event_loop.run(move |event, _, control_flow| match event {
        Event::DeviceEvent {
            event: DeviceEvent::MouseMotion { delta: (x, y) },
            ..
        } => state.mouse_movement = Vector2::new(x, y).map(|x| x as f32),
        Event::WindowEvent {
            window_id,
            ref event,
        } if window_id == window.id() => match event {
            WindowEvent::CloseRequested
            | WindowEvent::KeyboardInput {
                input:
                    KeyboardInput {
                        state: ElementState::Pressed,
                        virtual_keycode: Some(VirtualKeyCode::Escape),
                        ..
                    },
                ..
            } => {
                *control_flow = ControlFlow::Exit;
            }
            WindowEvent::Resized(new_size) => state.resize(*new_size),
            WindowEvent::ScaleFactorChanged { new_inner_size, .. } => {
                state.resize(**new_inner_size)
            }
            _ => state.input(event),
        },
        Event::RedrawRequested(window_id) if window_id == window.id() => {
            state.update();

            match state.render() {
                Ok(_) => {}
                Err(SurfaceError::Lost) => state.resize(state.size),
                Err(SurfaceError::OutOfMemory) => *control_flow = ControlFlow::Exit,
                Err(e) => error!("{:?}", e),
            }
        }
        Event::MainEventsCleared => {
            window.request_redraw();
        }
        _ => {}
    });
}
