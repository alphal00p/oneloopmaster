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
        .append_fun(&Slot::Out(0), "symbolica_sqrt", &[Slot::Param(0)], false)
        .unwrap();
    let app = translator.compile().unwrap().seal().unwrap();
    println!("fast_complex={packed}");

    let mut failures = 0;
    for (input, expected) in [
        (
            Complex::new(0.0, -4.0),
            Complex::new(2.0_f64.sqrt(), -2.0_f64.sqrt()),
        ),
        (Complex::new(-4.0, -0.0), Complex::new(0.0, -2.0)),
        (
            Complex::new(0.0, 4.0),
            Complex::new(2.0_f64.sqrt(), 2.0_f64.sqrt()),
        ),
        (Complex::new(-4.0, 0.0), Complex::new(0.0, 2.0)),
    ] {
        let actual = app.evaluate_single(&[input]);
        println!("sqrt({input:?}): actual={actual:?}, expected={expected:?}");
        let agrees = (actual - expected).norm() <= 1e-14;
        if !agrees {
            failures += 1;
        }
    }
    assert_eq!(
        failures, 0,
        "principal complex square-root branch is incorrect"
    );
}
