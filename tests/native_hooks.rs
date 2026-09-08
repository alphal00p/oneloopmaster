//! Bare Symbol callbacks, without transparent function-map definitions.
#[path = "support/fixtures.rs"]
mod fixtures;
use oneloop::{NativeEvaluator, ScalarIntegral};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, DoubleFloat, RealLike, SingleFloat},
};

#[test]
fn double_float_master_hooks_preserve_the_numeric_domain() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            oneloop::initialize().unwrap();
            let fixtures = fixtures::parse(include_str!("data/benchmark.txt"));
            for (family, symbol) in [
                (ScalarIntegral::A0, oneloop::A0()),
                (ScalarIntegral::B0, oneloop::B0()),
                (ScalarIntegral::DB0, oneloop::dB0()),
                (ScalarIntegral::C0, oneloop::C0()),
                (ScalarIntegral::D0, oneloop::D0()),
            ] {
                let parameters = (0..family.arity())
                    .map(|i| {
                        symbolica::symbol!(format!("native_hook_{}_{i}", family.name())).to_atom()
                    })
                    .collect::<Vec<_>>();
                let calls = [0, -1, -2].map(|tag| {
                    let mut arguments = vec![Atom::num(tag)];
                    arguments.extend(parameters.iter().cloned());
                    symbol.call(&arguments)
                });
                let prototype = DoubleFloat::from(0.0);
                let mut hooks = Atom::evaluator_multiple(&calls, &parameters)
                    .direct_translation(true)
                    .build()
                    .unwrap()
                    .map_coeff(&|z| {
                        Complex::new(
                            prototype.from_rational(&z.re),
                            prototype.from_rational(&z.im),
                        )
                    });
                let mut native = NativeEvaluator::<DoubleFloat>::new(family).unwrap();
                for row in fixtures.iter().filter(|r| r.family == family) {
                    let arguments = row
                        .args
                        .iter()
                        .map(|z| Complex::new(z.re.into(), z.im.into()))
                        .collect::<Vec<_>>();
                    let mut expected = [Complex::new(0.0.into(), 0.0.into()); 3];
                    let mut actual = expected;
                    native.evaluate(&arguments, &mut expected).unwrap();
                    hooks.evaluate(&arguments, &mut actual);
                    assert_eq!(actual, expected, "{} {}", family.name(), row.name);
                    row.check(&actual.map(|z| Complex::new(z.re.to_f64(), z.im.to_f64())));
                }
                if family == ScalarIntegral::A0 {
                    let precise = DoubleFloat::from(1.0) + DoubleFloat::from(2_f64.powi(-70));
                    let mut actual = [Complex::new(0.0.into(), 0.0.into()); 3];
                    hooks.evaluate(
                        &[
                            Complex::new(precise, 0.0.into()),
                            Complex::new(1.0.into(), 0.0.into()),
                        ],
                        &mut actual,
                    );
                    assert_eq!(actual[1].re, precise);
                    assert_ne!(
                        actual[1].re,
                        DoubleFloat::from(1.0),
                        "hook narrowed through f64"
                    );
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
