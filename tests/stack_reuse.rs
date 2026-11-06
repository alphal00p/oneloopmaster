//! `evaluator_with_stack_reuse` only reassigns the root's temporary slots: its
//! C0 and D0 values equal those of `evaluator` bit for bit, in binary64 and at
//! 128 bits, at the same exact inputs.
use oneloop::OneLoopExpressions;
use symbolica::{
    atom::Atom,
    domains::{
        float::{Complex, Float, FloatLike, SingleFloat},
        rational::Rational,
    },
};

fn parameters(prefix: &str, n: usize) -> Vec<Atom> {
    (0..n)
        .map(|i| symbolica::symbol!(format!("stack_reuse::{prefix}_{i}")).to_atom())
        .collect()
}

fn check(series: &[Atom], x: &[Atom], points: &[Vec<i64>]) {
    let e = OneLoopExpressions::new();
    let plain = e.evaluator(series, x).unwrap();
    let reused = e.evaluator_with_stack_reuse(series, x).unwrap();
    let to_f64 = |v: &Complex<Rational>| Complex::new(v.re.to_f64(), v.im.to_f64());
    let (mut a, mut b) = (plain.clone().map_coeff(&to_f64), reused.clone().map_coeff(&to_f64));
    let bits = 128;
    let p = Float::new(bits);
    let to_float =
        |v: &Complex<Rational>| Complex::new(p.from_rational(&v.re), p.from_rational(&v.im));
    let (mut fa, mut fb) = (
        plain.map_coeff_with_prec(&to_float, bits),
        reused.map_coeff_with_prec(&to_float, bits),
    );
    for point in points {
        let args: Vec<_> = point.iter().map(|&v| Complex::new(v as f64, 0.0)).collect();
        let (mut oa, mut ob) = (vec![Complex::new(0.0, 0.0); 3], vec![Complex::new(0.0, 0.0); 3]);
        a.evaluate(&args, &mut oa);
        b.evaluate(&args, &mut ob);
        for (x, y) in oa.iter().zip(&ob) {
            assert_eq!(x.re.to_bits(), y.re.to_bits(), "{point:?}: {oa:?} vs {ob:?}");
            assert_eq!(x.im.to_bits(), y.im.to_bits(), "{point:?}: {oa:?} vs {ob:?}");
        }
        let args: Vec<_> = point
            .iter()
            .map(|&v| Complex::new(Float::with_val(bits, v), Float::with_val(bits, 0)))
            .collect();
        let zero = Complex::new(p.zero(), p.zero());
        let (mut oa, mut ob) = (vec![zero.clone(); 3], vec![zero; 3]);
        fa.evaluate(&args, &mut oa);
        fb.evaluate(&args, &mut ob);
        assert_eq!(oa, ob, "{point:?}");
        for (x, y) in oa.iter().zip(&ob) {
            assert_eq!(x.re.get_precision(), y.re.get_precision(), "{point:?}");
            assert_eq!(x.im.get_precision(), y.im.get_precision(), "{point:?}");
        }
    }
}

#[test]
fn stack_reuse_keeps_c0_and_d0_values() {
    std::thread::Builder::new()
        .stack_size(512 << 20)
        .spawn(|| {
            let e = OneLoopExpressions::new();
            let x = parameters("c0", 7);
            let c0 = e.c0([&x[0], &x[1], &x[2]], [&x[3], &x[4], &x[5]], &x[6]);
            check(
                c0.coefficients(),
                &x,
                &[
                    vec![-1, -5, -4, 2, 3, 5, 1],
                    vec![-4000, 0, 0, 29929, 29929, 29929, 10000],
                    vec![3, 7, -2, 0, 1, 0, 2],
                ],
            );
            let x = parameters("d0", 11);
            let d0 = e.d0(
                [&x[0], &x[1], &x[2], &x[3], &x[4], &x[5]],
                [&x[6], &x[7], &x[8], &x[9]],
                &x[10],
            );
            check(
                d0.coefficients(),
                &x,
                &[
                    vec![0, 0, 0, 0, 7, -3, 0, 0, 0, 0, 1],
                    vec![1, 2, 3, 4, -10, -6, 5, 5, 5, 5, 2],
                ],
            );
        })
        .unwrap()
        .join()
        .unwrap();
}
