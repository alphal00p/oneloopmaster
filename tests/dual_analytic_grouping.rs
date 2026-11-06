//! Dual compilation replaces the exact grouping identity `if(z, z, 0)` by a
//! transparent identity, so a zero primal keeps every Taylor component, while a
//! genuine branch and the arithmetic grouping of the identity are retained.
//! Ported from iregi-rust's `crates/iregi-masters/tests/physical_derivatives.rs`.
use oneloop::OneLoopExpressions;
use symbolica::{
    atom::{Atom, Symbol},
    domains::{
        dual::HyperDual,
        float::{Complex, Float, FloatLike, SingleFloat},
        rational::Rational,
    },
};

const BITS: u32 = 384;
const SHAPE: [[usize; 2]; 6] = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]];
type C = Complex<Float>;

fn c(re: i64, im: i64) -> C {
    C::new(Float::with_val(BITS, re), Float::with_val(BITS, im))
}

fn map(v: &Complex<Rational>) -> C {
    C::new(
        v.re.to_multi_prec_float(BITS),
        v.im.to_multi_prec_float(BITS),
    )
}

#[test]
fn analytic_grouping_preserves_zero_primal_jets() {
    let context = OneLoopExpressions::new();
    let shape = || HyperDual::<f64>::new(SHAPE.map(Vec::from).to_vec());
    let x = symbolica::symbol!("dual_grouping_x").to_atom();
    let expression = Symbol::IF.call((&x, &x, 0));
    let mut evaluator = context
        .dual_evaluator(&[expression], std::slice::from_ref(&x), shape())
        .unwrap()
        .map_coeff_with_prec(&map, BITS)
        .into_batch_evaluator();
    let input = vec![c(0, 0), c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10)];
    let mut actual = vec![c(0, 0); SHAPE.len()];
    evaluator.evaluate_batch(&input, &mut actual, 1).unwrap();
    assert_eq!(
        actual, input,
        "an algebraic identity must preserve every Taylor component at zero"
    );

    // A genuine branch is not an analytic identity and must remain untouched.
    let conditional = Symbol::IF.call((&x, &x * Atom::num(2), 0));
    let mut conditional = context
        .dual_evaluator(&[conditional], std::slice::from_ref(&x), shape())
        .unwrap()
        .map_coeff_with_prec(&map, BITS)
        .into_batch_evaluator();
    conditional.evaluate_batch(&input, &mut actual, 1).unwrap();
    assert!(actual.iter().all(|v| v.re.is_zero() && v.im.is_zero()));

    // Retain the arithmetic grouping used to compare precisely equal large
    // masses while keeping an independently represented tiny invariant.
    let y = symbolica::symbol!("dual_grouping_y").to_atom();
    let p = symbolica::symbol!("dual_grouping_p").to_atom();
    let difference = &x - &y;
    let grouped = Symbol::IF.call((&difference, &difference, 0)) - &p;
    let mut grouped = context
        .dual_evaluator(&[grouped], &[x, y, p], shape())
        .unwrap()
        .map_coeff_with_prec(&map, BITS)
        .into_batch_evaluator();
    let mut arguments = vec![c(0, 0); 3 * SHAPE.len()];
    arguments[0] = c(29929, 0);
    arguments[SHAPE.len()] = c(29929, 0);
    arguments[1] = c(1, 0);
    arguments[SHAPE.len() + 2] = c(1, 0);
    let tiny = c(2, 0).pow(114).inv();
    arguments[2 * SHAPE.len()] = tiny.clone();
    grouped.evaluate_batch(&arguments, &mut actual, 1).unwrap();
    assert_eq!(actual[0], -tiny);
    assert_eq!(actual[1], c(1, 0));
    assert_eq!(actual[2], c(-1, 0));
    assert!(actual[3..].iter().all(|v| v.re.is_zero() && v.im.is_zero()));
}
