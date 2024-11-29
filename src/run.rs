use crate::*;

pub fn run(solver: &mut Solver, n_it: usize, n_out: usize) {
    println!("Starting Run");
    let mut iter = 0;

    println!("  writing first snapshot {:06}", iter);
    //solver.moments();
    solver.write_vtk(0);

    while iter < n_it {
        println!("  iter: {}", iter);
        let write_output = n_out > 0 && iter % n_out == 0;
        println!("    streaming...");
        solver.streaming();
        println!("    moments...");
        solver.moments();
        println!("    collision...");
        solver.collision();
        println!("    apply_bcs...");
        solver.apply_bcs();

        if write_output {
            println!("    writing snapshot {:06}", iter);
            solver.write_vtk(iter);
        }

        iter += 1;
    }
}
