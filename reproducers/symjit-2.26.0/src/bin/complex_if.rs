use symjit::{Complex, Composer, Config, Slot, Translator};

fn main() {
    let join_only = std::env::args().nth(1).as_deref() == Some("join-only");
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_fast_complex(true);
    config.set_direct(false);
    config.set_opt_level(2);
    config.set_simd(false);
    config.set_threads(false);
    config.set_simd_branch(false);

    let mut t = Translator::new(config);
    t.set_num_params(1);
    t.append_constant(Complex::new(7.0, 0.0)).unwrap();
    t.append_constant(Complex::new(11.0, 0.0)).unwrap();
    if join_only {
        t.append_join(
            &Slot::Out(0),
            &Slot::Param(0),
            &Slot::Const(0),
            &Slot::Const(1),
        )
        .unwrap();
    } else {
        t.append_if_else(&Slot::Param(0), 0).unwrap();
        t.append_assign(&Slot::Temp(0), &Slot::Const(0)).unwrap();
        t.append_goto(1).unwrap();
        t.append_label(0).unwrap();
        t.append_assign(&Slot::Temp(1), &Slot::Const(1)).unwrap();
        t.append_label(1).unwrap();
        t.append_join(
            &Slot::Out(0),
            &Slot::Param(0),
            &Slot::Temp(0),
            &Slot::Temp(1),
        )
        .unwrap();
    }
    let app = t.compile().unwrap().seal().unwrap();
    let mut failures = 0;
    for input in [
        Complex::new(0.0, 1.0),
        Complex::new(0.0, -1.0),
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
    ] {
        let expected = Complex::new(
            if input.re != 0.0 || input.im != 0.0 {
                7.0
            } else {
                11.0
            },
            0.0,
        );
        let actual = app.evaluate_single(&[input]);
        println!("if({input:?}, 7, 11): actual={actual:?}, expected={expected:?}");
        failures += usize::from(actual != expected);
    }
    assert_eq!(failures, 0, "complex truth must test both components");
}
