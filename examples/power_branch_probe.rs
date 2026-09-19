//! Diagnostic for exact-axis powers through native, JIT and portable evaluators.
use symbolica::evaluate::ExpressionEvaluator;
use symbolica::prelude::*;

fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(run)
        .unwrap()
        .join()
        .unwrap();
}

fn run() {
    let z = symbol!("power_probe_z").to_atom();
    let s = symbol!("power_probe_s").to_atom();
    let m0 = symbol!("power_probe_m0").to_atom();
    let m1 = symbol!("power_probe_m1").to_atom();
    let d = (&m1 - &m0 - &s).pow(2) - Atom::num(4) * &s * &m0;
    let half = Atom::num(1) / 2;
    let labels = [
        "sqrt(z)", "z^(-1/2)", "log(z)", "z^-1", "z^-2", "z^-3", "d", "sqrt(d)", "d^(-1/2)", "dB0",
    ];
    let expressions = [
        z.sqrt(),
        z.pow(-&half),
        z.log(),
        z.pow(-1),
        z.pow(-2),
        z.pow(-3),
        d.clone(),
        d.sqrt(),
        d.pow(-&half),
        oneloop::db0(&s, &m0, &m1, &Atom::num(1)).coefficients()[0].clone(),
    ];
    let exact = Atom::evaluator_multiple(&expressions, &[z, s, m0, m1])
        .direct_translation(true)
        .build()
        .unwrap();
    let mut reference = exact
        .clone()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let rows = [
        (-3., 0.),
        (-3., -0.),
        (3., 0.),
        (3., -0.),
        (-3., 4.),
        (-3., -4.),
        (3., 4.),
        (3., -4.),
        (0., 4.),
        (0., -4.),
        (-0., 4.),
        (-0., -4.),
    ]
    .map(|(re, im): (f64, f64)| {
        [
            Complex::new(re, im),
            Complex::new(3., 0.0_f64.copysign(im)),
            Complex::new(1., 0.),
            Complex::new(1., 0.),
        ]
    });
    let mut expected = vec![Complex::new(0., 0.); labels.len() * rows.len()];
    for (input, out) in rows.iter().zip(expected.chunks_exact_mut(labels.len())) {
        reference.evaluate(input, out);
        assert!((out[9].re - 0.47279971743743027).abs() < 1e-13);
    }
    for (mode, simd, fast_complex) in [
        ("scalar", false, false),
        ("SIMD", true, false),
        ("fast", false, true),
    ] {
        let settings = oneloop::jit_settings()
            .with_option("use_simd", simd.to_string())
            .with_option("fast_complex", fast_complex.to_string());
        let source = exact
            .clone()
            .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
        let mut original = source.jit_compile(settings.clone()).unwrap();
        let bytes = bincode::encode_to_vec(&source, bincode::config::standard()).unwrap();
        let (source, consumed): (ExpressionEvaluator<Complex<f64>>, usize) =
            bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
        assert_eq!(consumed, bytes.len());
        let mut restored = source.jit_compile(settings).unwrap();
        for (stage, compiled) in [("original", &mut original), ("restored", &mut restored)] {
            for size in [1, 3, 4, 5, 8, rows.len()] {
                let input = rows[..size].iter().flatten().copied().collect::<Vec<_>>();
                let mut actual = vec![Complex::new(f64::NAN, f64::NAN); size * labels.len()];
                compiled.batch_evaluate(&input, &mut actual, size);
                let mut expected = expected[..actual.len()].to_vec();
                // Equivalent arithmetic can produce opposite signed zeros in d.
                // Check each root against the native primitive at its actual base;
                // the complete dB0 must remain invariant under this root relabeling.
                for (actual, expected) in actual
                    .chunks_exact(labels.len())
                    .zip(expected.chunks_exact_mut(labels.len()))
                {
                    expected[7] = actual[6].sqrt();
                    expected[8] = actual[6].powf(&Complex::new(-0.5, 0.));
                }
                for (index, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
                    let error = *actual - *expected;
                    assert!(
                        error.re.hypot(error.im) < 1e-12,
                        "{mode}/{stage}/batch{size} input={:?} {}: {actual:?} != {expected:?}",
                        rows[index / labels.len()],
                        labels[index % labels.len()]
                    );
                }
                for out in actual.chunks_exact(labels.len()) {
                    for (root, inverse) in [(out[0], out[1]), (out[7], out[8])] {
                        let error = root * inverse - Complex::new(1., 0.);
                        assert!(error.re.hypot(error.im) < 1e-13);
                    }
                }
            }
        }
    }
}
