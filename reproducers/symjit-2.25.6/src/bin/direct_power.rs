use symjit::{Complex, Composer, Config, Slot, Translator};

fn main() {
    let direct = std::env::args().nth(1).as_deref() != Some("indirect");
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_direct(direct);
    config.set_opt_level(2);
    config.set_simd(false);
    config.set_threads(false);

    let mut t = Translator::new(config);
    t.set_num_params(1);
    t.append_constant(Complex::new(0.25, 0.0)).unwrap();
    t.append_powf(&Slot::Out(0), &Slot::Param(0), &Slot::Const(0), false)
        .unwrap();
    println!("compiling z^(1/4), direct={direct}");
    let app = t
        .compile()
        .expect("valid fractional power must compile")
        .seal()
        .unwrap();
    let actual = app.evaluate_single(&[Complex::new(16.0, 0.0)]);
    println!("actual={actual:?}, expected=2+0i");
    assert!((actual - Complex::new(2.0, 0.0)).norm() < 1e-14);
}
