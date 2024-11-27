use crate::*;
use lattice::*;

// Coord Iterator
pub fn coord_iter(aabb: AABB<3>) -> impl std::iter::Iterator<Item = Coord<3>> {
    let size = box_buffer_size(&aabb);
    (0..size)
        .map(move |index| 
            linear_to_coord_in_box(index, &aabb)
        )
}

pub struct Solver {
    grid_dimensions: AABB<3>,
    distributions: Array4D,
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
    pub fn new(grid_dimensions: AABB<3>, omega: f32, c_sqr: f32, inflow_density: f32, inflow_accel: f32) -> Self {
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
        result.init();
        result
    }

    fn init(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            for q in 0..27 {
                let value = self.inflow_density * D3Q27_W[q];
                self.distributions.set_q(&coord, q as i32, value)
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
                let mut pressure = 0.0 ;
                let mut u = Vec3::zero(); 
                for q_i in 0..27 {
                    let q = self.distributions.get_q(&coord, q_i);
                    pressure += q;
                    u += self.directions[q_i as usize] * q;
                }
                u /= pressure;
                self.pressure.set(&coord, pressure);
                self.velocity.set(&coord, u);
        }

    }

    pub fn collision(&mut self) {
        for coord in coord_iter(self.grid_dimensions) {
            let u = self.velocity.get(&coord);
            let p = self.pressure.get(&coord);
            for q_i in 0..27 {
                // Calculate equilibrium
                let dir = self.directions[q_i as usize];
                let dir_u = dir.dot(&u);
                let w_i = D3Q27_W[q_i as usize]; 

                let t1 = (3.0 * dir_u) / self.c_sqr;
                let t2 = (9.0 * dir_u * dir_u) / (2.0 * self.c_sqr * self.c_sqr);
                let t3 = - (3.0 * u.dot(&u)) / (2.0 * self.c_sqr);
                let q_eq = w_i * p * (1.0 + t1 + t2 + t3); 
                
                // relax
                let q = self.distributions.get_q(&coord, q_i);
                let new_q = q + self.omega * (q_eq - q);
                self.distributions.set_q(&coord, q_i, new_q);
            }
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

    pub fn apply_inflow_bc(&mut self, coord &Coord<3>) {
        for q_i in 0..27 {
            let t = self.inflow_density * self.inflow_accel * D3Q27_W[q_i];
        }
    }

    pub fn apply_bcs(&mut self) { 

        // Inflow,
        let mut inflow_face = self.grid_dimensions;
        inflow_face[(2, 1)] = inflow_face[(2, 0)];

        // Outlfow
        let mut outflow_face = self.grid_dimensions;
        outflow_face[(2, 0)] = inflow_face[(2, 1)];

        // Bounce back
        let mut bb_0 = self.grid_dimensions;
        bb_0[(0, 1)] = bb_0[(0, 0)];
        for coord in coord_iter(bb_0) {
            self.apply_bounce_back(&coord);
        }

        let mut bb_1 = self.grid_dimensions;
        bb_1[(0, 0)] = bb_0[(0, 1)];
        for coord in coord_iter(bb_1) {
            self.apply_bounce_back(&coord);
        }

        let mut bb_2 = self.grid_dimensions;
        bb_2[(1, 1)] = bb_2[(1, 0)];
        for coord in coord_iter(bb_2) {
            self.apply_bounce_back(&coord);
        }

        let mut bb_3 = self.grid_dimensions;
        bb_3[(1, 0)] = bb_3[(1, 1)];
        for coord in coord_iter(bb_3) {
            self.apply_bounce_back(&coord);
        }
    }

    pub fn write_vtk(&self, _i: usize) {

    }
}
