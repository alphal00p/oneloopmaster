use symjit::{Complex, Composer, Config, Slot, Translator};

fn main() {
    let packed = std::env::args().nth(1).as_deref() != Some("generic");
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_fast_complex(packed);
    config.set_direct(false);
    config.set_opt_level(2);
    config.set_simd(false);
    config.set_threads(false);

    let mut translator = Translator::new(config);
    translator.set_num_params(1);
    translator
        .append_pow(&Slot::Out(0), &Slot::Param(0), -1, false)
        .unwrap();
    let app = translator.compile().unwrap().seal().unwrap();
    println!("fast_complex={packed}");

    let mut failures = 0;
    for (input, expected) in [
        (Complex::new(2.0, 0.0), Complex::new(0.5, 0.0)),
        (Complex::new(1e300, 0.0), Complex::new(1e-300, 0.0)),
        (Complex::new(0.0, -1e300), Complex::new(0.0, 1e-300)),
        (Complex::new(1e300, 1e300), Complex::new(5e-301, -5e-301)),
        (Complex::new(1e-300, 0.0), Complex::new(1e300, 0.0)),
        (Complex::new(0.0, -1e-300), Complex::new(0.0, 1e300)),
        (Complex::new(1e-300, 1e-300), Complex::new(5e299, -5e299)),
    ] {
        let actual = app.evaluate_single(&[input]);
        println!("reciprocal({input:?}): actual={actual:?}, expected={expected:?}");
        let agrees = [(actual.re, expected.re), (actual.im, expected.im)]
            .into_iter()
            .all(|(a, e)| {
                a.is_finite()
                    && if e == 0.0 {
                        a == 0.0
                    } else {
                        (a / e - 1.0).abs() <= 1e-14
                    }
            });
        if !agrees {
            failures += 1;
        }
    }
    assert_eq!(
        failures, 0,
        "complex reciprocal lost a representable component"
    );
}
