//! The smaller build must preserve automatic evaluation and reject unavailable backends.
#[path = "support/fixtures.rs"]
mod fixtures;
use oneloop::{EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral};
use symbolica::domains::float::{Complex, Float};

#[test]
fn default_backend_preserves_scalar_and_batch_fixtures() {
    let expected = if cfg!(feature = "generated-evaluators") {
        EvaluationBackend::Native
    } else {
        EvaluationBackend::Expression
    };
    assert_eq!(oneloop::DEFAULT_BACKEND, expected);
    let fixtures = fixtures::parse(include_str!("data/parity.txt"));
    assert_eq!(fixtures.len(), 293);
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        let rows = fixtures
            .iter()
            .filter(|row| row.family == family)
            .collect::<Vec<_>>();
        for row in &rows {
            let mut output = [Complex::new(0., 0.); 3];
            oneloop::evaluate(family, &row.args, &mut output).unwrap();
            row.check(&output);
        }
        let input = rows
            .iter()
            .flat_map(|row| row.args.iter().copied())
            .collect::<Vec<_>>();
        let mut output = vec![Complex::new(f64::NAN, f64::NAN); rows.len() * 3];
        oneloop::evaluate_batch(family, &input, &mut output, rows.len()).unwrap();
        for (row, values) in rows.iter().zip(output.chunks_exact(3)) {
            row.check(values);
        }
        let mut invalid = input.clone();
        invalid.last_mut().unwrap().re = -1.;
        let sentinel = vec![Complex::new(17., -9.); output.len()];
        output.clone_from(&sentinel);
        assert!(oneloop::evaluate_batch(family, &invalid, &mut output, rows.len()).is_err());
        assert_eq!(output, sentinel);
    }
}

#[test]
fn explicit_native_selection_matches_feature_availability() {
    assert_eq!(
        EvaluationBackend::Native.is_available(),
        cfg!(feature = "generated-evaluators")
    );
    assert!(EvaluationBackend::Expression.is_available());
    assert!(EvaluationBackend::SymJit.is_available());
    let family = ScalarIntegral::A0;
    let native = NativeEvaluator::<f64>::new(family);
    let native_float = NativeEvaluator::<Float>::with_binary_precision(family, 128);
    let precise = PrecisionEvaluator::with_binary_precision_and_backend(
        family,
        128,
        EvaluationBackend::Native,
    );
    let mut output = [Complex::new(17., -9.); 3];
    let manual = oneloop::evaluate_with_backend(
        family,
        &[Complex::new(2., 0.), Complex::new(1., 0.)],
        &mut output,
        EvaluationBackend::Native,
    );
    for result in [
        native.map(|_| ()),
        native_float.map(|_| ()),
        precise.map(|_| ()),
        manual,
    ] {
        if cfg!(feature = "generated-evaluators") {
            result.unwrap();
        } else {
            assert!(result.unwrap_err().contains("generated-evaluators"));
        }
    }
    if !cfg!(feature = "generated-evaluators") {
        assert_eq!(output, [Complex::new(17., -9.); 3]);
    }
    assert_eq!(
        PrecisionEvaluator::with_binary_precision(family, 128)
            .unwrap()
            .backend(),
        oneloop::DEFAULT_BACKEND
    );
}
