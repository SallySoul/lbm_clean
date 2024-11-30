use crate::*;
use lattice::*;
use nalgebra::{matrix, vector};
use vtkio::model::*;
use rand::prelude::*;
use rand::distributions::{Distribution, Uniform};



// Coord Iterator
pub fn coord_iter(aabb: AABB<3>) -> impl std::iter::Iterator<Item = Coord<3>> {
    let size = box_buffer_size(&aabb);
    (0..size).map(move |index| linear_to_coord_in_box(index, &aabb))
}

pub fn cell_coord_iter(aabb: AABB<3>) -> impl std::iter::Iterator<Item = Coord<3>> {
    let mut cell_bounds = aabb.clone();
    cell_bounds.set_column(1, &cell_bounds.column(1).add_scalar(-1));
    let size = box_buffer_size(&cell_bounds);
    (0..size).map(move |index| linear_to_coord_in_box(index, &cell_bounds))
}

pub fn cell_count(aabb: AABB<3>) -> usize {
    let mut cell_bounds = aabb.clone();
    cell_bounds.set_column(1, &cell_bounds.column(1).add_scalar(-1));
    box_buffer_size(&cell_bounds)
}

pub struct Solver {
    grid_dimensions: AABB<3>,
    pub distributions: Array4D,
    distributions_buffer: Array4D,
    pressure: Array3D,
    velocity: VelArray,
    offsets: [Coord<3>; 27],
    directions: [Vec3; 27],
    omega: f32,
    c_sqr: f32,
    inflow_density: f32,
    inflow_accel: f32,
}

impl Solver {
    pub fn new(
        grid_dimensions: AABB<3>,
        omega: f32,
        c_sqr: f32,
        inflow_density: f32,
        inflow_accel: f32,
    ) -> Self {
        let q_bounds = nalgebra::matrix![0, 26];
        #[allow(clippy::toplevel_ref_arg)]
        let dimensions = nalgebra::stack![grid_dimensions; q_bounds];

        let mut result = Solver {
            grid_dimensions,
            distributions: Array4D::new(dimensions),
            distributions_buffer: Array4D::new(dimensions),
            pressure: Array3D::new(grid_dimensions),
            velocity: VelArray::new(grid_dimensions),
            offsets: gen_d3q27_offsets(),
            directions: gen_d3q27_directions(),
            omega,
            c_sqr,
            inflow_density,
            inflow_accel,
        };
        result
    }

    pub fn equilibrium_init(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            for q in 0..27 {
                let value = self.inflow_density * D3Q27_W[q];
                self.distributions.set_q(&coord, q as i32, value)
            }
        }
    }

    pub fn flow_init(&mut self) {

        let mut rng = rand::thread_rng();
        let dist = Uniform::from(0.0..0.1);
        for coord in coord_iter(self.grid_dimensions) {
            for q_i in 0..27 {
                let value = self.inflow_density * dist.sample(&mut rng);
                self.distributions.set_q(&coord, q_i, value);
            }

            if coord[1] > 20 {
                let q_i = 6;
                let w_i = D3Q27_W[q_i];
                let w = 1.0 / w_i;
                let value = self.inflow_density * w;
                self.distributions.set_q(&coord, q_i as i32, value);
            } else {
                let q_i = 5;
                let w_i = D3Q27_W[q_i];
                let w = 1.0 / w_i;
                let value = self.inflow_density * w;
                self.distributions.set_q(&coord, q_i as i32, value);
            }
        }
    }

    pub fn streaming(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            for q_i in 0..27 {
                // Get q value
                let q = self.distributions.get_q(&coord, q_i);

                // Get neighbor intex
                let neighbor_coord = coord + self.offsets[q_i as usize];

                if box_contains_coord(&self.grid_dimensions, &neighbor_coord) {
                    self.distributions_buffer.set_q(&neighbor_coord, q_i, q);
                }
            }
        }
        std::mem::swap(&mut self.distributions, &mut self.distributions_buffer);
    }

    pub fn moments(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            let mut pressure = 0.0;
            let mut u = Vec3::zero();
            for q_i in 0..27 {
                let q = self.distributions.get_q(&coord, q_i);
                pressure += q;
                u += self.directions[q_i as usize] * q;
            }
            if pressure.abs() > 0.00001 {
                u /= pressure;
            } 
            self.pressure.set(&coord, pressure);
            self.velocity.set(&coord, u);
        }
    }

    pub fn collision(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            let u = self.velocity.get(&coord);
            let p = self.pressure.get(&coord);
            let mut i_d = 0.0;
            let mut n_d = 0.0;
            let mut eq_d = 0.0;
            for q_i in 0..27 {
                // Calculate equilibrium
                let dir = self.directions[q_i as usize];
                let dir_u = dir.dot(&u);
                let w_i = D3Q27_W[q_i as usize];

                let t1 = (3.0 * dir_u) / self.c_sqr;
                let t2 = (9.0 * dir_u * dir_u) / (2.0 * self.c_sqr * self.c_sqr);
                let t3 = -(3.0 * u.dot(&u)) / (2.0 * self.c_sqr);
                let q_eq = w_i * p * (1.0 + t1 + t2 + t3);

                // relax
                let q = self.distributions.get_q(&coord, q_i);
                let new_q = q + self.omega * (q_eq - q);
                self.distributions.set_q(&coord, q_i, new_q);

                i_d += q;
                n_d += new_q;
                eq_d += q_eq;
            }
            println!("q: {}, eq_q: {}, n_d: {}, p: {}, u: {:?}", i_d, eq_d, n_d, p, u);
        }
    }

    pub fn apply_bounce_back(&mut self, coord: &Coord<3>) {
        let mut new_q = [0.0; 27];
        for q_i in 0..27 {
            let q = self.distributions.get_q(coord, q_i);
            new_q[D3Q27_OPP[q_i as usize]] = q;
        }
        for q_i in 0..27 {
            self.distributions.set_q(coord, q_i, new_q[q_i as usize]);
        }
    }
/*
    pub fn apply_inflow_bc(&mut self, coord &Coord<3>) {
        for q_i in 0..27 {
            let t = self.inflow_density * self.inflow_accel * D3Q27_W[q_i];
        }


    }

    pub fn apply_outflow_bc(&mut self, coord: &Coord<3>) {
        for q_i in 0..27 {
            let t = self.inflow_density * self.inflow_accel * D3Q27_W[q_i];
  k      }
    }
*/
    pub fn apply_bcs(&mut self) {
        let x_min = self.grid_dimensions[(0, 0)];
        let x_max = self.grid_dimensions[(0, 1)];
        let y_min = self.grid_dimensions[(1, 0)];
        let y_max = self.grid_dimensions[(1, 1)];
        let z_min = self.grid_dimensions[(2, 0)];
        let z_max = self.grid_dimensions[(2, 1)];

        // Bottom
        let bb_0 = matrix![x_min, x_max; y_min, y_max; z_min, z_max;];
        for coord in coord_iter(bb_0) {
            self.apply_bounce_back(&coord);
        }

        // Top
        let bb_1 = matrix![x_min, x_max; y_max, y_max; z_min, z_max;];
        for coord in coord_iter(bb_1) {
            self.apply_bounce_back(&coord);
        }

        // Left (x_min)
        let bb_2 = matrix![x_min, x_min; y_min + 1, y_max - 1; z_min, z_max;];
        for coord in coord_iter(bb_2) {
            self.apply_bounce_back(&coord);
        }

        // Right
        let bb_3 = matrix![x_max, x_max; y_min + 1, y_max - 1; z_min, z_max;];
        for coord in coord_iter(bb_3) {
            self.apply_bounce_back(&coord);
        }

        // front
        let bb_4 = matrix![x_min + 1, x_max - 1; y_min + 1, y_max - 1;  z_min, z_min];
        for coord in coord_iter(bb_4) {
            self.apply_bounce_back(&coord);
        }

        let bb_5 = matrix![x_min + 1, x_max - 1; y_min + 1, y_max - 1; z_max, z_max];
        for coord in coord_iter(bb_5) {
            self.apply_bounce_back(&coord);
        }

        /*
        // Inflow,
        let mut inflow_face = self.grid_dimensions;
        inflow_face[(2, 1)] = inflow_face[(2, 0)];

        // Outlfow
        let mut outflow_face = self.grid_dimensions;
        outflow_face[(2, 0)] = inflow_face[(2, 1)];

        */
        /*
        // Bounce back
        let mut bb_0 = self.grid_dimensions;
        bb_0[(0, 1)] = bb_0[(0, 0)];
        for coord in coord_iter(bb_0) {
            println!("bb_coord: {:?}", coord);
            self.apply_bounce_back(&coord);
        }

        let mut bb_1 = self.grid_dimensions;
        bb_1[(0, 0)] = bb_0[(0, 1)];
        for coord in coord_iter(bb_1) {
            println!("bb_coord: {:?}", coord);
            self.apply_bounce_back(&coord);
        }

        let mut bb_2 = self.grid_dimensions;
        bb_2[(1, 1)] = bb_2[(1, 0)];
        for coord in coord_iter(bb_2) {
            println!("bb_coord: {:?}", coord);
            self.apply_bounce_back(&coord);
        }

        let mut bb_3 = self.grid_dimensions;
        bb_3[(1, 0)] = bb_3[(1, 1)];
        for coord in coord_iter(bb_3) {
            println!("bb_coord: {:?}", coord);
            self.apply_bounce_back(&coord);
        }
        */
    }

    pub fn write_vtk(&self, i: usize) {
        let buffer_size = box_buffer_size(&self.grid_dimensions);
        //let distributions = vec![vec![0.0; buffer_size]; 27];
        let mut density = Vec::with_capacity(buffer_size);
        let mut velocity = Vec::with_capacity(3 * buffer_size);
        let mut points = Vec::with_capacity(3 * buffer_size);
        let mut qs = vec![Vec::with_capacity(buffer_size); 27];
        for coord in coord_iter(self.grid_dimensions) {
            points.push(coord[0] as f32);
            points.push(coord[1] as f32);
            points.push(coord[2] as f32);

            density.push(self.pressure.get(&coord));

            let vel = self.velocity.get(&coord);
            velocity.push(vel[0]);
            velocity.push(vel[1]);
            velocity.push(vel[2]);

            for q_i in 0..27 {
                let q = self.distributions.get_q(&coord, q_i);
                qs[q_i as usize].push(q);
            }
        }

        let n_cells = cell_count(self.grid_dimensions);
        let mut connectivity = Vec::with_capacity(n_cells);
        let mut offsets = Vec::with_capacity(n_cells);
        let mut cell_types = Vec::with_capacity(n_cells);
        let mut offset = 8;
        for cell_coord in cell_coord_iter(self.grid_dimensions) {
            let n_1 = cell_coord + vector![0, 0, 1];
            let n_2 = cell_coord + vector![0, 1, 0];
            let n_3 = cell_coord + vector![1, 0, 0];
            let n_4 = cell_coord + vector![0, 1, 1];
            let n_5 = cell_coord + vector![1, 1, 0];
            let n_6 = cell_coord + vector![1, 0, 1];
            let n_7 = cell_coord + vector![1, 1, 1];

            let vertices = [&cell_coord, &n_3, &n_6, &n_1, &n_2, &n_5, &n_7, &n_4];
            for v in vertices {
                let index = coord_to_linear_in_box(v, &self.grid_dimensions) as u64;
                connectivity.push(index);
            }

            offsets.push(offset);
            cell_types.push(CellType::Hexahedron);
            offset += 8;
        }

        let mut point_attributes = vec![
            Attribute::DataArray(DataArrayBase {
                name: format!("density"),
                elem: ElementType::Scalars {
                    num_comp: 1,
                    lookup_table: None,
                },
                data: IOBuffer::F32(density),
            }),
            Attribute::DataArray(DataArrayBase {
                name: format!("velocity"),
                elem: ElementType::Scalars {
                    num_comp: 3,
                    lookup_table: None,
                },
                data: IOBuffer::F32(velocity),
            }),
        ];

        for q_i in 0..27 {
            point_attributes.push(Attribute::DataArray(DataArrayBase {
                name: format!("q_{}", q_i),
                elem: ElementType::Scalars {
                    num_comp: 1,
                    lookup_table: None,
                },
                data: IOBuffer::F32(qs[q_i].clone()),
            }));
        }

        Vtk {
            version: Version { major: 1, minor: 0 },
            title: String::new(),
            byte_order: ByteOrder::LittleEndian,
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                points: IOBuffer::F32(points),
                cells: Cells {
                    cell_verts: VertexNumbers::XML {
                        connectivity,
                        offsets,
                    },
                    types: cell_types,
                },
                data: Attributes {
                    point: point_attributes,
                    cell: vec![],
                },
            }),
        }
        .export(format!("vtk_test/data_{:06}.vtu", i))
        .unwrap();
    }
}
