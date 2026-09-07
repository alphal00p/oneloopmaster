//! Run with `cargo run --example basic` after the README dependency bootstrap.
use symbolica::prelude::*;

fn main() {
    // Construct every Symbolica value inside this thread. This also accommodates
    // larger C0/D0 examples without activating a second Symbolica worker.
    std::thread::Builder::new()
        .name("oneloop-example".into())
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let p = parse!("p");
            let mass_squared = Atom::num(1);
            let mu_squared = Atom::num(1);
            let tadpole = oneloop::a0(&mass_squared, &mu_squared);
            println!("A0(m²=1, mu²=1), exact Laurent coefficients:");
            for (power, coefficient) in [0, -1, -2].into_iter().zip(tadpole.coefficients()) {
                println!("  epsilon^{power}: {coefficient}");
            }

            let bubble = oneloop::b0(&p, &mass_squared, &mass_squared, &mu_squared);
            let views = bubble.coefficients().each_ref().map(Atom::as_view);
            let mut evaluator = Atom::evaluator_multiple(&views, &[p])
                .build()
                .expect("compile B0")
                .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
            let mut output = [Complex::new(0.0, 0.0); 3];
            evaluator.evaluate(&[Complex::new(-1.0, 0.0)], &mut output);
            println!("B0(p²=-1, m0²=m1²=mu²=1):");
            println!("  [finite, simple pole, double pole] = {output:?}");
            let expected = [-0.1520447048200202, 1.0, 0.0];
            for (value, reference) in output.into_iter().zip(expected) {
                assert!((value.re - reference).hypot(value.im) < 1e-12);
            }
        })
        .expect("start example thread")
        .join()
        .expect("example thread failed");
}
