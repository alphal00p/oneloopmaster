use symjit::{Complex, Composer, Config, Slot, Translator};

fn main() {
    let packed = std::env::args().nth(1).as_deref() == Some("packed");
    let mut config = Config::default();
    config.set_complex(true);
    config.set_fastmath(false);
    config.set_fast_complex(packed);
    config.set_direct(false);
    config.set_opt_level(2);
    config.set_simd(false);
    config.set_threads(false);

    let mut t = Translator::new(config);
    t.set_num_params(1);
    t.append_fun(&Slot::Temp(0), "symbolica_conj", &[Slot::Param(0)], false)
        .unwrap();
    t.append_mul(&Slot::Out(0), &[Slot::Param(0), Slot::Temp(0)], 0)
        .unwrap();
    let app = t.compile().unwrap().seal().unwrap();
    let z = Complex::new(0.1, 0.2);
    let actual = app.evaluate_single(&[z]);
    // Separate binary64 products cancel exactly in the imaginary component.
    let expected = z * z.conj();
    println!("fastmath=false, fast_complex={packed}");
    println!("z={z:?}; z*conj(z): actual={actual:?}, expected={expected:?}");
    assert_eq!(
        actual.im, 0.0,
        "strict complex multiplication introduced an imaginary residue"
    );
}
