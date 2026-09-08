//! Expected callback errors must not poison a shared native family cache.
use oneloop::{NativeEvaluator, NativeFloat, ScalarIntegral};
use std::panic::{AssertUnwindSafe, catch_unwind};
use symbolica::{
    atom::Atom,
    domains::float::{Complex, DoubleFloat, Float},
};

fn check_domain<T: NativeFloat>(bits: u32, convert: impl Fn(f64) -> T) {
    let value = |real, imaginary| Complex::new(convert(real), convert(imaginary));
    for (family, symbol, momenta) in [
        (ScalarIntegral::A0, oneloop::A0(), 0),
        (ScalarIntegral::B0, oneloop::B0(), 1),
        (ScalarIntegral::DB0, oneloop::dB0(), 1),
        (ScalarIntegral::C0, oneloop::C0(), 3),
        (ScalarIntegral::D0, oneloop::D0(), 6),
    ] {
        let info = symbol.get_evaluation_info().expect("master EvaluationInfo");
        let callbacks = [0, -1, -2].map(|tag| {
            info.get_evaluator::<Complex<T>>(&[Atom::num(tag).as_view()])
                .expect("native numeric callback")
        });
        let mut valid = vec![value(1., 0.); family.arity()];
        valid[..momenta].fill(value(0., 0.));
        let mut expected = core::array::from_fn::<_, 3, _>(|_| value(0., 0.));
        NativeEvaluator::<T>::with_binary_precision(family, bits)
            .unwrap()
            .evaluate(&valid, &mut expected)
            .unwrap();

        let mut upper_mass = valid.clone();
        upper_mass[momenta].im = convert(0.125);
        let mut zero_scale = valid.clone();
        *zero_scale.last_mut().unwrap() = value(0., 0.);
        let mut negative_scale = valid.clone();
        *negative_scale.last_mut().unwrap() = value(-1., 0.);

        // Keep a partially consumed valid group before the first rejection.
        assert_eq!(callbacks[0](&valid), expected[0]);
        for invalid in [upper_mass, zero_scale, negative_scale] {
            for (tag, callback) in callbacks.iter().enumerate() {
                let panic = catch_unwind(AssertUnwindSafe(|| callback(&invalid)))
                    .expect_err("invalid master domain must panic");
                let message = panic
                    .downcast_ref::<String>()
                    .map(String::as_str)
                    .or_else(|| panic.downcast_ref::<&str>().copied())
                    .expect("descriptive callback panic");
                assert!(
                    message.contains("valid master arguments") && !message.contains("poison"),
                    "{} tag {tag}, {bits} bits: {message}",
                    family.name()
                );
                // Exercise the same family, numeric domain and precision
                // immediately after catching the error. For binary64, the
                // subsequent invalid tags also detect accidentally cached
                // placeholder outputs from the rejected point.
                for (valid_tag, callback) in callbacks.iter().enumerate() {
                    assert_eq!(
                        callback(&valid),
                        expected[valid_tag],
                        "{} tag {valid_tag}, {bits} bits after rejection",
                        family.name()
                    );
                }
            }
        }
    }
}

#[test]
fn invalid_native_hook_domains_leave_all_numeric_caches_usable() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            oneloop::initialize().unwrap();
            check_domain::<f64>(53, |value| value);
            check_domain::<DoubleFloat>(106, DoubleFloat::from);
            for bits in [128, 256] {
                check_domain::<Float>(bits, |value| Float::with_val(bits, value));
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
