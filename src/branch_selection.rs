//! Select analytic regions without substituting the parametric formulas.
use std::convert::Infallible;

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::float::{Complex, Float, SingleFloat},
    id::{Pattern, ReplaceWith, Replacement},
    transformer::Transformer,
    utils::BorrowedOrOwned,
};

/// Select native Symbolica `if(condition, nonzero, zero)` branches at a sample point.
///
/// The native replacements are applied, in their supplied order, to a temporary
/// copy of each condition only. Neither the retained branch nor an undecidable
/// condition is numerically substituted. Nested conditionals are treated at every
/// depth, including inside conditions, function arguments, bases and exponents.
/// Unselected branches are not traversed.
///
/// Symbolica performs both pattern matching and numerical predicate evaluation.
/// Integer, rational and complex-rational literal conditions are decided exactly;
/// other evaluable conditions use its arbitrary-precision complex arithmetic, without a zero tolerance or
/// conversion to machine floats. Unknown and non-finite conditions remain `if`s.
/// Floating-point replacement values are first represented by their exact binary
/// rational values, avoiding premature rounding during predicate substitution.
/// This does not recover precision already lost when constructing a value or
/// inside a caller's replacement callback.
/// Selecting a region does not make the resulting formula valid in other regions.
///
/// ```
/// use oneloop::select_branch;
/// use symbolica::prelude::*;
/// let x = symbol!("select_branch_example_x").to_atom();
/// let expression = Symbol::IF.call((&x, &x + 1, &x - 1));
/// assert_eq!(select_branch(&expression, &[Replacement::new(x.clone(), 2)]), x + 1);
/// ```
pub fn select_branch(expression: impl AtomCore, replacement_rules: &[Replacement]) -> Atom {
    let rules = replacement_rules
        .iter()
        .map(exact_replacement)
        .collect::<Vec<_>>();
    select_branch_impl(
        expression,
        |condition| Ok::<_, Infallible>(condition.replace_multiple(&rules)),
        true,
    )
    .unwrap_or_else(|never| match never {})
}

fn exact_float_coefficients(expression: AtomView<'_>) -> Atom {
    expression.map_coefficient(|coefficient| {
        if let CoefficientView::Float(real, imaginary) = coefficient
            && let (Some(real), Some(imaginary)) = (
                real.to_float().try_to_rational(),
                imaginary.to_float().try_to_rational(),
            )
        {
            Coefficient::Complex(Complex::new(real, imaginary))
        } else {
            coefficient.to_owned()
        }
    })
}

fn protected_rhs(expression: AtomView<'_>) -> Atom {
    let exact = exact_float_coefficients(expression);
    // This native, numerically evaluable identity keeps substituted numbers
    // out of symbolic arithmetic until predicate evaluation. In particular,
    // forming sqrt(1 + an exact tiny float) here would otherwise ask Symbolica
    // to factor a potentially thousands-of-bits integer during normalization.
    // PI is a built-in nonzero constant, and IF evaluates its active arm lazily.
    // The guard lives only in the temporary sample, never in the result.
    Symbol::IF.call((Symbol::PI, &exact, &exact))
}

pub(crate) fn exact_replacement(replacement: &Replacement) -> Replacement {
    let mut replacement = replacement.clone();
    replacement.rhs = match replacement.rhs {
        ReplaceWith::Pattern(pattern) => {
            let mut pattern = pattern.yield_owned();
            exact_pattern_values(&mut pattern);
            ReplaceWith::Pattern(BorrowedOrOwned::Owned(pattern))
        }
        ReplaceWith::Map(map) => ReplaceWith::Map(Box::new(move |matches| {
            protected_rhs(map(matches).as_view())
        })),
    };
    replacement
}

fn exact_pattern_values(pattern: &mut Pattern) {
    let mut work = vec![pattern];
    while let Some(pattern) = work.pop() {
        match pattern {
            Pattern::Literal(value) => *value = protected_rhs(value.as_view()),
            Pattern::Wildcard(..) => {}
            Pattern::Fn(_, arguments)
            | Pattern::Mul(arguments)
            | Pattern::Add(arguments)
            | Pattern::Alternative(arguments) => work.extend(arguments.iter_mut()),
            Pattern::Pow(arguments) => work.extend(arguments.iter_mut()),
            Pattern::Transformer(transformer) => {
                if let Some(input) = &mut transformer.0 {
                    work.push(input);
                }
                // Preserve the caller's full native transformer. Its output is
                // protected before the surrounding condition is normalized.
                transformer
                    .1
                    .push(Transformer::Map(Box::new(|value, _, out| {
                        *out = protected_rhs(value);
                        Ok(())
                    })));
            }
        }
    }
}

/// Variant of [`select_branch`] for host-language native replacement adapters.
///
/// The callback receives only a condition, and its result is used only to decide
/// that condition. Callback errors are propagated; numerical inability to decide a condition is not an
/// error and leaves the original parametric condition in place.
/// Unlike [`select_branch`], this low-level adapter cannot prepare replacement
/// right-hand sides: the callback must itself avoid premature sample rounding.
pub fn select_branch_with<E>(
    expression: impl AtomCore,
    replace_condition: impl FnMut(AtomView<'_>) -> Result<Atom, E>,
) -> Result<Atom, E> {
    select_branch_impl(expression, replace_condition, false)
}

fn select_branch_impl<E>(
    expression: impl AtomCore,
    mut replace_condition: impl FnMut(AtomView<'_>) -> Result<Atom, E>,
    cache: bool,
) -> Result<Atom, E> {
    // An explicit work stack avoids a helper-imposed recursion limit, and lets
    // us decide each IF before entering either arm (native IF is lazy).
    let mut work = vec![Work::Visit(expression.as_atom_view())];
    let mut values = Vec::new();
    let mut memo = std::collections::HashMap::new();
    let mut predicates = std::collections::HashMap::new();
    while let Some(next) = work.pop() {
        match next {
            Work::Visit(view) => {
                if cache && !matches!(view, AtomView::Num(_) | AtomView::Var(_)) {
                    if let Some(value) = memo.get(&view) {
                        values.push(Atom::clone(value));
                        continue;
                    }
                    work.push(Work::Save(view));
                }
                match view {
                    AtomView::Num(_) | AtomView::Var(_) => values.push(view.to_owned()),
                    AtomView::Fun(fun)
                        if fun.get_symbol() == Symbol::IF && fun.get_nargs() == 3 =>
                    {
                        let mut arguments = fun.iter();
                        let condition = arguments.next().unwrap();
                        let yes = arguments.next().unwrap();
                        let no = arguments.next().unwrap();
                        // Probe the original predicate before stripping grouping
                        // IFs inside it. Rebuilding x-|x| after removing if(x,x,0)
                        // can change floating-point evaluation order at exact zero.
                        if cache {
                            let key = condition.to_owned();
                            let truth = if let Some(truth) = predicates.get(&key) {
                                *truth
                            } else {
                                let sample = replace_condition(condition)?;
                                let truth = numeric_truth(sample.as_view());
                                predicates.insert(key, truth);
                                truth
                            };
                            if let Some(truth) = truth {
                                work.push(Work::Visit(if truth { yes } else { no }));
                                continue;
                            }
                        }
                        work.push(Work::Condition { yes, no });
                        work.push(Work::Visit(condition));
                    }
                    AtomView::Fun(fun) => {
                        schedule(&mut work, Build::Function(fun.get_symbol()), fun.iter());
                    }
                    AtomView::Pow(power) => {
                        let (base, exponent) = power.get_base_exp();
                        schedule(&mut work, Build::Power, [base, exponent].into_iter());
                    }
                    AtomView::Mul(product) => schedule(&mut work, Build::Product, product.iter()),
                    AtomView::Add(sum) => schedule(&mut work, Build::Sum, sum.iter()),
                }
            }
            Work::Save(view) => {
                memo.insert(view, values.last().unwrap().clone());
            }
            Work::Condition { yes, no } => {
                let condition = values.pop().unwrap();
                let truth = if let Some(truth) = predicates.get(&condition) {
                    *truth
                } else {
                    let sample = replace_condition(condition.as_view())?;
                    let truth = numeric_truth(sample.as_view());
                    if cache {
                        predicates.insert(condition.clone(), truth);
                    }
                    truth
                };
                match truth {
                    Some(true) => work.push(Work::Visit(yes)),
                    Some(false) => work.push(Work::Visit(no)),
                    None => {
                        // Keep the parametric condition, not the substituted
                        // sample, while still selecting IFs inside both arms.
                        values.push(condition);
                        work.push(Work::Build(Build::Function(Symbol::IF), 3));
                        work.push(Work::Visit(no));
                        work.push(Work::Visit(yes));
                    }
                }
            }
            Work::Build(kind, count) => {
                let mut children = values.split_off(values.len() - count);
                values.push(match kind {
                    Build::Function(symbol) => symbol.call(&children),
                    Build::Power => {
                        let exponent = children.pop().unwrap();
                        children.pop().unwrap().pow(exponent)
                    }
                    Build::Product => children.into_iter().product(),
                    Build::Sum => children.into_iter().sum(),
                });
            }
        }
    }
    Ok(values.pop().unwrap())
}

enum Work<'a> {
    Visit(AtomView<'a>),
    Save(AtomView<'a>),
    Condition { yes: AtomView<'a>, no: AtomView<'a> },
    Build(Build, usize),
}

enum Build {
    Function(Symbol),
    Power,
    Product,
    Sum,
}

fn schedule<'a>(
    work: &mut Vec<Work<'a>>,
    kind: Build,
    children: impl Iterator<Item = AtomView<'a>>,
) {
    let children = children.collect::<Vec<_>>();
    work.push(Work::Build(kind, children.len()));
    work.extend(children.into_iter().rev().map(Work::Visit));
}

pub(crate) fn numeric_truth(condition: AtomView<'_>) -> Option<bool> {
    if let AtomView::Num(number) = condition {
        return match number.get_coeff_view() {
            CoefficientView::Natural(..) | CoefficientView::Large(..) => Some(!number.is_zero()),
            CoefficientView::Float(real, imaginary)
                if real.to_float().is_finite() && imaginary.to_float().is_finite() =>
            {
                Some(!number.is_zero())
            }
            _ => None,
        };
    }

    let bits = numeric_precision(condition)?;
    let value: Complex<Float> = condition
        .evaluate_with_prec::<Atom, _>(&std::collections::HashMap::new(), bits)
        .ok()?;
    value.is_finite().then(|| !value.is_zero())
}

pub(crate) fn numeric_precision(condition: AtomView<'_>) -> Option<u32> {
    // Keep all supplied coefficient precision. Extra working bits reduce
    // cancellation in numerical predicates; they do not constitute a proof at
    // a singular boundary, nor recover digits absent from an input float.
    let mut bits = 256u32;
    let mut pending = vec![condition];
    while let Some(view) = pending.pop() {
        match view {
            AtomView::Num(number) => match number.get_coeff_view().to_owned() {
                Coefficient::Float(value) => {
                    bits = bits.max(value.re.prec().max(value.im.prec()).saturating_add(64));
                }
                Coefficient::Complex(value) => {
                    for part in [&value.re, &value.im] {
                        let size = part
                            .numerator_ref()
                            .significant_bits()
                            .max(part.denominator_ref().significant_bits());
                        bits = bits.max(u32::try_from(size).ok()?.saturating_add(64));
                    }
                }
                // These may live in an inactive native IF arm. Let the lazy
                // evaluator reject them only if that arm is actually needed.
                _ => {}
            },
            AtomView::Fun(fun) => pending.extend(fun.iter()),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                pending.extend([base, exponent]);
            }
            AtomView::Mul(product) => pending.extend(product.iter()),
            AtomView::Add(sum) => pending.extend(sum.iter()),
            AtomView::Var(_) => {}
        }
    }
    Some(bits)
}
