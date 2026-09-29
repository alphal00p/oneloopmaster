use symjit::{Complex, Composer, Config, Defuns, Slot, Translator};

fn callback(args: &[Complex<f64>]) -> Complex<f64> {
    args[0] * Complex::new(2., -3.) + args[1] * Complex::new(-5., 7.)
}

fn main() {
    let simd = std::env::args().nth(1).as_deref() != Some("scalar");
    let mut defuns = Defuns::new();
    defuns
        .add_sliced_func("affine", Box::new(callback))
        .unwrap();
    let mut config = Config::default();
    config.set_defuns(defuns);
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_fast_complex(false);
    config.set_direct(true);
    config.set_opt_level(2);
    config.set_simd(simd);
    config.enable_simd512(false);
    config.set_threads(false);
    let mut translator = Translator::new(config);
    translator.set_num_params(2);
    translator
        .append_external_fun(&Slot::Out(0), "affine", &[Slot::Param(0), Slot::Param(1)])
        .unwrap();
    let app = translator.compile().unwrap().seal().unwrap();
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
                Complex::new(1. + i as f64, 11. + i as f64),
                Complex::new(21. + i as f64, 31. + i as f64),
            ]
        })
        .collect();
    let expected: Vec<_> = inputs.chunks_exact(2).map(callback).collect();
    for (row, expected) in inputs.chunks_exact(2).zip(&expected) {
        assert_eq!(app.evaluate_single(row), *expected);
    }
    let mut actual = vec![Complex::new(f64::NAN, f64::NAN); 4];
    app.evaluate_matrix(&inputs, &mut actual, 4);
    println!("actual={actual:?}\nexpected={expected:?}");
    assert_eq!(
        actual, expected,
        "scalar complex callback lanes were repacked incorrectly"
    );
}
