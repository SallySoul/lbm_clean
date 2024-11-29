use crate::*;

/// x, y, z, q
pub struct Array4D {
    dimensions: AABB<4>,
    size: usize,
    pub buffer: Vec<f32>,
}

impl Array4D {
    pub fn new(dimensions: AABB<4>) -> Self {
        let size = box_buffer_size(&dimensions);

        Array4D {
            dimensions,
            size,
            buffer: vec![0.0; size],
        }
    }

    pub fn get(&self, coord: &Coord<4>) -> f32 {
        let index = coord_to_linear_in_box(coord, &self.dimensions);
        self.buffer[index]
    }

    pub fn get_q(&self, grid_coord: &Coord<3>, q: i32) -> f32 {
        #[allow(clippy::toplevel_ref_arg)]
        let coord = nalgebra::stack![grid_coord; nalgebra::vector![q]];
        let index = coord_to_linear_in_box(&coord, &self.dimensions);
        self.buffer[index]
    }

    pub fn set(&mut self, coord: Coord<4>, value: f32) {
        let index = coord_to_linear_in_box(&coord, &self.dimensions);
        self.buffer[index] = value;
    }

    pub fn set_q(&mut self, grid_coord: &Coord<3>, q: i32, value: f32) {
        #[allow(clippy::toplevel_ref_arg)]
        let coord = nalgebra::stack![grid_coord; nalgebra::vector![q]];
        let index = coord_to_linear_in_box(&coord, &self.dimensions);
        self.buffer[index] = value;
    }
}

pub struct Array3D {
    dimensions: AABB<3>,
    size: usize,
    buffer: Vec<f32>,
}

impl Array3D {
    pub fn new(dimensions: AABB<3>) -> Self {
        let size = box_buffer_size(&dimensions);

        Array3D {
            dimensions,
            size,
            buffer: vec![0.0; size],
        }
    }

    pub fn get(&self, coord: &Coord<3>) -> f32 {
        let index = coord_to_linear_in_box(coord, &self.dimensions);
        self.buffer[index]
    }

    pub fn set(&mut self, coord: &Coord<3>, value: f32) {
        let index = coord_to_linear_in_box(&coord, &self.dimensions);
        self.buffer[index] = value;
    }
}

pub struct VelArray {
    dimensions: AABB<3>,
    size: usize,
    buffer: Vec<Vec3>,
}

impl VelArray {
    pub fn new(dimensions: AABB<3>) -> Self {
        let size = box_buffer_size(&dimensions);

        VelArray {
            dimensions,
            size,
            buffer: vec![Vec3::zero(); size],
        }
    }

    pub fn get(&self, coord: &Coord<3>) -> Vec3 {
        let index = coord_to_linear_in_box(coord, &self.dimensions);
        self.buffer[index]
    }

    pub fn set(&mut self, coord: &Coord<3>, value: Vec3) {
        let index = coord_to_linear_in_box(&coord, &self.dimensions);
        self.buffer[index] = value;
    }
}
