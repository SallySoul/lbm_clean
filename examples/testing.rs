use lbm_clean::*;
use nalgebra::matrix;

fn main() {
    let grid_dimensions = matrix![0, 30; 0, 40; 0, 50];
    //let omega = 1.85;
    let omega = 0.03;
    let c_sqr = 1.0 / 3.0;
    let inflow_density = 0.1;
    let inflow_accel = 0.015;
    let mut solver = Solver::new(grid_dimensions, omega, c_sqr, inflow_density, inflow_accel);
    run(&mut solver, 40, 1);
}
