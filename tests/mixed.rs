use oneloop::OneLoopExpressions;
use symbolica::{atom::Atom, parse, prelude::Complex};

#[test]
fn triangle_and_box_share_one_function_map() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let context = OneLoopExpressions::new();
            let x = parse!("mixed_scale");
            let p = [
                -&x,
                -Atom::num(2) * &x,
                -Atom::num(3) * &x,
                -Atom::num(4) * &x,
                -Atom::num(12) * &x,
                -Atom::num(15) * &x,
            ];
            let m = [parse!("1/2"), parse!("7/10"), parse!("7/5"), parse!("2")];
            let c = context.c0([&p[0], &p[1], &p[2]], [&m[0], &m[1], &m[2]], &Atom::num(1));
            let d = context.d0(
                [&p[0], &p[1], &p[2], &p[3], &p[4], &p[5]],
                [&m[0], &m[1], &m[2], &m[3]],
                &Atom::num(1),
            );
            let combined =
                Atom::num(2) * &c.coefficients()[0] + Atom::num(3) * &d.coefficients()[0];
            let mut evaluator = context
                .evaluator(
                    &[
                        c.coefficients()[0].clone(),
                        d.coefficients()[0].clone(),
                        combined,
                    ],
                    &[x],
                )
                .unwrap()
                .map_coeff(&|v| Complex::new(v.re.to_f64(), v.im.to_f64()));
            for scale in [0.5, 1., 2.] {
                let mut out = [Complex::new(0., 0.); 3];
                evaluator.evaluate(&[Complex::new(scale, 0.)], &mut out);
                assert!(out.iter().all(|z| z.re.is_finite() && z.im.is_finite()));
                let difference = out[2] - out[0] * 2. - out[1] * 3.;
                assert!(difference.re.hypot(difference.im) < 1e-12);
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
