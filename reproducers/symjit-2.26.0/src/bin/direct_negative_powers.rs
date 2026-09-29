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
    t.append_pow(&Slot::Out(0), &Slot::Param(0), -2, false)
        .unwrap();
    t.append_pow(&Slot::Out(1), &Slot::Param(0), -3, false)
        .unwrap();
    let app = t.compile().unwrap().seal().unwrap();
    let mut actual = [Complex::new(0.0, 0.0); 2];
    app.evaluate(&[Complex::new(2.0, 0.0)], &mut actual);
    let expected = [Complex::new(0.25, 0.0), Complex::new(0.125, 0.0)];
    println!("direct={direct}: [2^(-2), 2^(-3)]={actual:?}, expected={expected:?}");
    assert_eq!(actual, expected);
}
