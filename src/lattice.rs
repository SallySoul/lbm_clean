use crate::*;
use nalgebra::vector;

pub fn gen_d3q27_directions() -> [Vec3; 27] {
    [
        vector![0.0, 0.0, 0.0],
        vector![1.0, 0.0, 0.0],
        vector![-1.0, 0.0, 0.0],
        vector![0.0, 1.0, 0.0],
        vector![0.0, -1.0, 0.0],
        vector![0.0, 0.0, 1.0],
        vector![0.0, 0.0, -1.0],
        vector![1.0, 1.0, 0.0],
        vector![1.0, -1.0, 0.0],
        vector![-1.0, 1.0, 0.0],
        vector![-1.0, -1.0, 0.0],
        vector![1.0, 0.0, 1.0],
        vector![1.0, 0.0, -1.0],
        vector![-1.0, 0.0, 1.0],
        vector![-1.0, 0.0, -1.0],
        vector![0.0, 1.0, 1.0],
        vector![0.0, 1.0, -1.0],
        vector![0.0, -1.0, 1.0],
        vector![0.0, -1.0, -1.0],
        vector![1.0, 1.0, 1.0],
        vector![1.0, 1.0, -1.0],
        vector![1.0, -1.0, 1.0],
        vector![-1.0, 1.0, 1.0],
        vector![1.0, -1.0, -1.0],
        vector![-1.0, -1.0, 1.0],
        vector![-1.0, 1.0, -1.0],
        vector![-1.0, -1.0, -1.0],
    ]
}

pub fn gen_d3q27_offsets() -> [Coord<3>; 27] {
    [
        vector![0, 0, 0],   // 0
        vector![1, 0, 0],   // 1
        vector![-1, 0, 0],  // 2
        vector![0, 1, 0],   // 3
        vector![0, -1, 0],  // 4
        vector![0, 0, 1],   // 5
        vector![0, 0, -1],  // 6
        vector![1, 1, 0],   // 7
        vector![1, -1, 0],  // 8
        vector![-1, 1, 0],  // 9
        vector![-1, -1, 0], // 10
        vector![1, 0, 1],   // 11
        vector![1, 0, -1],  // 12
        vector![-1, 0, 1],  // 13
        vector![-1, 0, -1], // 14
        vector![0, 1, 1],   // 15
        vector![0, 1, -1],  // 16
        vector![0, -1, 1],  // 17
        vector![0, -1, -1], // 18
        vector![1, 1, 1],   // 19
        vector![1, 1, -1],  // 20
        vector![1, -1, 1],  // 21
        vector![-1, 1, 1],  // 22
        vector![1, -1, -1], // 23
        vector![-1, -1, 1], // 24
        vector![-1, 1, -1], // 25
        vector![-1, -1, -1], // 26
    ]
}

pub static D3Q27_OPP: [usize; 27] = [
    0,
    2,
    1,
    4,
    3,
    6,
    5,
    10,
    9,
    8,
    7,
    14,
    13,
    12,
    11,
    18,
    17,
    16,
    15,
    26,
    24,
    25,
    23,
    22,
    20,
    21,
    19,
];

pub static D3Q27_W: [f32; 27] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];

#[cfg(test)]
mod unit_tests {
    use super::*;
    
    #[test]
    fn opposites() {
        let offsets = gen_d3q27_offsets();
        for i in 0..27 {
            let o = D3Q27_OPP[i];
            let r: Coord<3> = offsets[i] + offsets[o];
            for d in 0..3 {
                assert_eq!(r[d], 0);
            }
        }
    }

    #[test]
    fn weights() {
        let s:f32 = D3Q27_W.iter().sum();
        assert!((1.0 - s).abs() < 0.00001);
    }

    #[test]
    fn dirs() {
        let mut sum = Vec3::zero();
        for d in gen_d3q27_directions() {
            sum += d;
        }
        for d in 0..3 {
            assert!(sum[d].abs() < 0.000001);
        }
    }
}
