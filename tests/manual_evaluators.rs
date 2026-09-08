//! Public manually composed JIT API, independent of the scalar-family wrappers.
use oneloop::{JitEvaluator, OneLoopExpressions};
use symbolica::prelude::*;

#[test]
fn manual_cache_clones_batches_and_constant_only_expressions() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let context = OneLoopExpressions::new();
            let x = symbol!("manual_cache_x").to_atom();
            let expressions = [2 * &x + 3, x.conj(), x.pow(-2)];
            let original = context.jit_evaluator(&expressions, &[x]).unwrap();
            assert_eq!((original.input_count(), original.output_count()), (1, 3));
            let bytes = original.to_bytes().unwrap();
            let header_end = bytes.iter().position(|&b| b == 0).unwrap() + 1;
            for offset in [header_end, header_end + 8] {
                let mut corrupt = bytes.clone();
                corrupt[offset..offset + 8].copy_from_slice(&0u64.to_le_bytes());
                assert!(JitEvaluator::from_bytes(&corrupt).is_err());
            }
            let mut cloned = original.clone();
            drop(original);
            let mut restored = JitEvaluator::from_bytes(&bytes).unwrap();
            let input = [(1., 0.), (2., -0.1), (-1., 0.3), (0., 2.), (3., -4.)]
                .map(|(re, im)| Complex::new(re, im));
            for evaluator in [&mut cloned, &mut restored] {
                evaluator.evaluate_batch(&[], &mut [], 0).unwrap();
                assert!(
                    evaluator
                        .evaluate(&[], &mut [Complex::new(0., 0.); 3])
                        .is_err()
                );
                assert!(evaluator.evaluate_batch(&[], &mut [], usize::MAX).is_err());
                for count in [1, 4, 5] {
                    let mut output = vec![Complex::new(f64::NAN, f64::NAN); 3 * count];
                    evaluator
                        .evaluate_batch(&input[..count], &mut output, count)
                        .unwrap();
                    for (&z, actual) in input[..count].iter().zip(output.chunks_exact(3)) {
                        let square = z * z;
                        let norm = square.re * square.re + square.im * square.im;
                        let expected = [
                            Complex::new(2. * z.re + 3., 2. * z.im),
                            Complex::new(z.re, -z.im),
                            Complex::new(square.re / norm, -square.im / norm),
                        ];
                        for (a, b) in actual.iter().zip(expected) {
                            assert!((a.re - b.re).hypot(a.im - b.im) < 1e-12, "{a:?} != {b:?}");
                        }
                    }
                }
            }
            let constant = context
                .jit_evaluator(&[Atom::num(2), Symbol::PI.to_atom()], &[])
                .unwrap();
            let mut constant = JitEvaluator::from_bytes(&constant.to_bytes().unwrap()).unwrap();
            assert_eq!((constant.input_count(), constant.output_count()), (0, 2));
            let mut output = [Complex::new(f64::NAN, f64::NAN); 10];
            constant.evaluate_batch(&[], &mut output, 5).unwrap();
            for row in output.chunks_exact(2) {
                assert_eq!(row[0], Complex::new(2., 0.));
                assert!((row[1].re - std::f64::consts::PI).abs() < 1e-15 && row[1].im == 0.);
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
