//! Public numeric callbacks for operations the upstream JIT lowers incorrectly.
use std::sync::OnceLock;
use symbolica::{
    atom::{Atom, AtomView, EvaluationInfo, Symbol},
    domains::float::{Complex, Real, SingleFloat},
};

pub(super) fn register() -> (Symbol, Symbol) {
    // Enter state before locking: OneLOop's state initializer registers these too.
    let _ = symbolica::get_symbol!("oneloopmaster::__jit_sqrt");
    static SYMBOLS: OnceLock<(Symbol, Symbol)> = OnceLock::new();
    *SYMBOLS.get_or_init(|| {
        let sqrt = symbolica::symbol!(
            "oneloopmaster::__jit_sqrt",
            eval = EvaluationInfo::new().register(|args: &[Complex<f64>]| args[0].sqrt())
        );
        let nonzero = symbolica::symbol!(
            "oneloopmaster::__jit_nonzero",
            eval = EvaluationInfo::new().register(|args: &[Complex<f64>]| {
                Complex::new(if args[0].is_zero() { 0. } else { 1. }, 0.)
            })
        );
        (sqrt, nonzero)
    })
}

pub(super) fn lower(view: AtomView<'_>, symbols: (Symbol, Symbol)) -> Option<Atom> {
    match view {
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            (exponent == Atom::num((1, 2))).then(|| symbols.0.call((base,)))
        }
        AtomView::Fun(call) if call.get_symbol() == Symbol::SQRT && call.get_nargs() == 1 => Some(
            symbols
                .0
                .call(&call.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>()),
        ),
        AtomView::Fun(call) if call.get_symbol() == Symbol::IF && call.get_nargs() == 3 => {
            let mut args = call.iter();
            Some(Symbol::IF.call((
                symbols.1.call((args.next().unwrap(),)),
                args.next().unwrap(),
                args.next().unwrap(),
            )))
        }
        _ => None,
    }
}
