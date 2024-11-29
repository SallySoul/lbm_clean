use lbm_clean::*;
use nalgebra::matrix;

fn main() {
    let grid_dimensions = matrix![0, 6; 0, 6; 0, 6];
    //let omega = 1.85;
    let omega = 0.8;
    let c_sqr = 1.0 / 3.0;
    let inflow_density = 0.1;
    let inflow_accel = 0.015;
    let mut solver = Solver::new(grid_dimensions, omega, c_sqr, inflow_density, inflow_accel);
    run(&mut solver, 20, 1);
}
