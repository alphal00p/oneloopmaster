//! Diagnostic of native complex branch conditions before and after portable JIT restoration.
use symbolica::evaluate::JITCompiledEvaluator;
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
    let z = symbol!("complex_branch_probe_z").to_atom();
    let labels = ["if(z)", "log", "conj", "if(z-conj(z))", "A0"];
    let expressions = [
        Symbol::IF.call((&z, 7, 11)),
        z.log(),
        z.conj(),
        Symbol::IF.call((&z - z.conj(), 7, 11)),
        oneloop::a0(&z, &Atom::num(1)).coefficients()[0].clone(),
    ];
    let exact = Atom::evaluator_multiple(&expressions, &[z])
        .direct_translation(true)
        .build()
        .unwrap();
    let mut reference = exact
        .clone()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let inputs = [
        Complex::new(0., 1.),
        Complex::new(0., -1.),
        Complex::new(2., -0.1),
        Complex::new(2., 0.),
        Complex::new(-2., 0.),
    ];
    let mut failures = 0;
    for (name, simd, branch, fast_complex) in [
        ("scalar", false, false, false),
        ("simd", true, false, false),
        ("simd_branch", true, true, false),
        ("fast_complex", false, false, true),
    ] {
        let settings = oneloop::jit_settings()
            .with_option("use_simd", simd.to_string())
            .with_option("simd_branch", branch.to_string())
            .with_option("fast_complex", fast_complex.to_string());
        let mut original = exact.jit_compile::<Complex<f64>>(settings.clone()).unwrap();
        let bytes = original.export_portable().unwrap();
        let mut restored =
            JITCompiledEvaluator::<Complex<f64>>::import_portable(&bytes, settings).unwrap();
        for (stage, compiled) in [("original", &mut original), ("restored", &mut restored)] {
            let mut expected = vec![Complex::new(0., 0.); inputs.len() * labels.len()];
            let mut actual = expected.clone();
            for (input, output) in inputs.iter().zip(expected.chunks_exact_mut(labels.len())) {
                reference.evaluate(&[*input], output);
            }
            for (input, output) in inputs.iter().zip(actual.chunks_exact_mut(labels.len())) {
                compiled.evaluate(&[*input], output);
            }
            for (layout, output) in [
                ("point", actual.clone()),
                ("batch", {
                    compiled.batch_evaluate(&inputs, &mut actual, inputs.len());
                    actual
                }),
            ] {
                for (index, (actual, expected)) in output.iter().zip(&expected).enumerate() {
                    let error = *actual - *expected;
                    let error_norm = error.re.hypot(error.im);
                    if !error_norm.is_finite() || error_norm >= 1e-12 {
                        eprintln!(
                            "{name}/{stage}/{layout} input={:?} {} actual={actual:?} expected={expected:?}",
                            inputs[index / labels.len()],
                            labels[index % labels.len()]
                        );
                        failures += 1;
                    }
                }
            }
        }
    }
    assert_eq!(failures, 0, "complex branch probe mismatches");
}
