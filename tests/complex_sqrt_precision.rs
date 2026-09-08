//! Regression for Numerica's exact half-angle constant precision provenance.
use symbolica::domains::float::{Complex, Float, Real};

#[test]
fn nearly_imaginary_sqrt_does_not_inherit_cancelled_real_precision() {
    for (bits, tiny, tolerance) in [(512, "1e-154", "1e-140"), (3456, "1e-1100", "1e-1000")] {
        let n = |text| Float::parse(text, Some(bits)).unwrap();
        let absolute = |value: Float| if value.is_negative() { -value } else { value };
        let expected = (n("0.5")).sqrt();
        for real_negative in [false, true] {
            for imaginary_negative in [false, true] {
                // Only two relative bits remain in the negligible real part;
                // its absolute uncertainty is far below the requested digits.
                let mut re = Float::parse(tiny, Some(2)).unwrap();
                let mut im = n("1");
                if real_negative {
                    re = -re;
                }
                if imaginary_negative {
                    im = -im;
                }
                let actual = Complex::new(re, im).sqrt();
                let expected_im = if imaginary_negative {
                    -expected.clone()
                } else {
                    expected.clone()
                };
                assert!(absolute(actual.re.clone() - expected.clone()) < n(tolerance));
                assert!(absolute(actual.im.clone() - expected_im) < n(tolerance));
                assert!(
                    actual.re.prec() >= bits - 8 && actual.im.prec() >= bits - 8,
                    "phase precision was poisoned: real={}, imaginary={}, working={bits}",
                    actual.re.prec(),
                    actual.im.prec()
                );
            }
        }
    }
}
