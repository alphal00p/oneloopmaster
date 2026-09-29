use symjit::{Applet, Complex, Composer, Config, Defuns, Slot, Translator};

fn config(simd: bool) -> Config {
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_fast_complex(false);
    config.set_direct(true);
    config.set_opt_level(2);
    config.set_simd(simd);
    config.set_simd_branch(false);
    config.enable_simd512(false);
    config.set_threads(false);
    config
}

fn conditional_child(simd: bool) -> Applet {
    let mut t = Translator::new(config(simd));
    t.set_num_params(2);
    t.append_if_else(&Slot::Param(0), 1).unwrap();
    t.append_assign(&Slot::Out(0), &Slot::Param(1)).unwrap();
    t.append_goto(2).unwrap();
    t.append_label(1).unwrap();
    t.append_add(&Slot::Out(0), &[Slot::Param(1), Slot::Param(1)], 0)
        .unwrap();
    t.append_label(2).unwrap();
    t.compile().unwrap().seal().unwrap()
}

fn main() {
    let simd = std::env::args().nth(1).as_deref() != Some("scalar");
    let mut defuns = Defuns::new();
    defuns.add_applet("child", conditional_child(simd));
    let mut settings = config(simd);
    settings.set_defuns(defuns);
    let mut t = Translator::new(settings);
    t.set_num_params(2);
    t.append_external_fun(&Slot::Temp(0), "child", &[Slot::Param(0), Slot::Param(1)])
        .unwrap();
    t.append_add(&Slot::Out(0), &[Slot::Temp(0), Slot::Param(1)], 0)
        .unwrap();
    let app = t.compile().unwrap().seal().unwrap();
    println!(
        "SIMD requested={simd}, compiled={}",
        app.compiled_simd.is_some()
    );
    if simd {
        assert!(app.compiled_simd.is_some(), "requires a SIMD-capable host");
    }
    let inputs: Vec<_> = (0..4)
        .flat_map(|i| {
            [
                Complex::new((i % 2) as f64, 0.),
                Complex::new(i as f64 + 0.25, i as f64 - 0.75),
            ]
        })
        .collect();
    let expected: Vec<_> = inputs
        .chunks_exact(2)
        .map(|row| row[1] * if row[0].re == 0. { 3. } else { 2. })
        .collect();
    for (row, expected) in inputs.chunks_exact(2).zip(&expected) {
        assert_eq!(app.evaluate_single(row), *expected);
    }
    let mut actual = vec![Complex::new(f64::NAN, f64::NAN); 4];
    app.evaluate_matrix(&inputs, &mut actual, 4);
    println!("actual={actual:?}\nexpected={expected:?}");
    assert_eq!(
        actual, expected,
        "nested SIMD bailout did not trigger scalar retry"
    );
}
