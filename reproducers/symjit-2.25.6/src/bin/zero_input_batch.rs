use symjit::{Complex, Composer, Config, Slot, Translator};

fn main() {
    let simd = std::env::args().nth(1).as_deref() != Some("scalar");
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_direct(false);
    config.set_opt_level(2);
    config.set_simd(simd);
    config.enable_simd512(false);
    config.set_simd_branch(false);
    config.set_threads(false);

    let mut t = Translator::new(config);
    t.set_num_params(0);
    t.append_constant(Complex::new(2.0, 0.0)).unwrap();
    t.append_constant(Complex::new(std::f64::consts::PI, 0.0))
        .unwrap();
    t.append_assign(&Slot::Out(0), &Slot::Const(0)).unwrap();
    t.append_assign(&Slot::Out(1), &Slot::Const(1)).unwrap();
    let app = t.compile().unwrap().seal().unwrap();
    let expected = [
        Complex::new(2.0, 0.0),
        Complex::new(std::f64::consts::PI, 0.0),
    ];
    let mut single = [Complex::new(0.0, 0.0); 2];
    app.evaluate(&[], &mut single);
    assert_eq!(single, expected);
    let mut output = [Complex::new(0.0, 0.0); 10];
    eprintln!(
        "calling evaluate_matrix with zero parameters, 2 outputs, 5 rows; SIMD requested={simd}, compiled={}",
        app.compiled_simd.is_some()
    );
    app.evaluate_matrix(&[], &mut output, 5);
    println!("actual={output:?}");
    for row in output.chunks_exact(2) {
        assert_eq!(row, expected);
    }
}
