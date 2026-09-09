//! Complete, inspectable master expressions, independent of numerical backends.
use crate::{LaurentSeries, ScalarIntegral};
use std::collections::{HashMap, HashSet};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    evaluate::FunctionMap,
    transcendental::TranscendentalFunctions,
};

/// Resource limits for fully substituting the transparent native definitions.
///
/// Generic boxes can be much larger than their compact FunctionMap form. Limits
/// fail explicitly; they never return a truncated formula or opaque helper.
#[derive(Clone, Copy, Debug)]
pub struct ExpressionOptions {
    /// Cumulative expression nodes visited or materialized, including memoized
    /// result copies. This is a work bound, not just the final expression size.
    pub max_nodes: usize,
    /// Maximum recursive expression/definition depth.
    pub max_depth: usize,
}

impl Default for ExpressionOptions {
    fn default() -> Self {
        Self {
            max_nodes: 1_000_000,
            max_depth: 512,
        }
    }
}

/// Return all complete master coefficients in `[epsilon^0, epsilon^-1,
/// epsilon^-2]` order from a master call, with `mu_squared` last.
///
/// For example, `get_expression(crate::B0().call((s, m0, m1, mu_squared)))`.
/// This inspection call omits the Laurent tag of a numerical coefficient call.
/// Composite arguments are kept as expressions; no conversion to a variable or
/// independent assumption declarations are required.
///
/// Every registered OneLOop helper is substituted by its exact native Symbolica
/// body. Input symbols retain their existing attributes: `Real` and `Positive`
/// can simplify conjugation and branch predicates, but no symbol is retagged.
/// Arbitrary string tags do not imply mathematical assumptions. Remaining `if`
/// expressions preserve branches that cannot be decided symbolically.
///
/// This is an inspection/explicit-compilation API, not the compact numerical
/// evaluation path. Inlining can duplicate arithmetic that the mapped evaluator
/// evaluates once; it does not inherit that evaluator's numerical conditioning.
pub fn get_expression(master: impl AtomCore) -> Result<LaurentSeries, String> {
    get_expression_with_options(master, ExpressionOptions::default())
}

/// Version of [`get_expression`] with explicit expansion resource limits.
pub fn get_expression_with_options(
    master: impl AtomCore,
    options: ExpressionOptions,
) -> Result<LaurentSeries, String> {
    let (family, arguments) = master_arguments(master)?;
    get_expression_for_family_with_options(family, &arguments, options)
}

/// Decode an inspection call such as `oneloopmaster::B0(s, m0, m1, mu_squared)`.
///
/// Arguments may be arbitrary Symbolica expressions. Inspection calls omit the
/// leading Laurent tag used by the numerical coefficient symbols: the result
/// contains all three coefficients. The `oneloop::` namespace is accepted as an
/// inspection shorthand; it is not a second set of numerical evaluation hooks.
pub fn master_arguments(master: impl AtomCore) -> Result<(ScalarIntegral, Vec<Atom>), String> {
    let AtomView::Fun(call) = master.as_atom_view() else {
        return Err("expected a scalar master call with its arguments, for example oneloopmaster::B0(psq,m0_squared,m1_squared,mu_squared)".into());
    };
    let family = match call.get_symbol().get_name() {
        "oneloopmaster::A0" | "oneloop::A0" => ScalarIntegral::A0,
        "oneloopmaster::B0" | "oneloop::B0" => ScalarIntegral::B0,
        "oneloopmaster::dB0" | "oneloop::dB0" => ScalarIntegral::DB0,
        "oneloopmaster::C0" | "oneloop::C0" => ScalarIntegral::C0,
        "oneloopmaster::D0" | "oneloop::D0" => ScalarIntegral::D0,
        _ => {
            return Err(format!(
                "{} is not a supported scalar master symbol",
                call.get_symbol()
            ));
        }
    };
    if call.get_nargs() != family.arity() {
        return Err(format!(
            "{} inspection expects {} arguments including mu_squared last and no Laurent tag; received {}",
            family.name(),
            family.arity(),
            call.get_nargs()
        ));
    }
    Ok((
        family,
        call.iter().map(|argument| argument.to_owned()).collect(),
    ))
}

/// Enum-based alternative to [`get_expression`]; all arguments are expressions.
pub fn get_expression_for_family(
    family: ScalarIntegral,
    arguments: &[Atom],
) -> Result<LaurentSeries, String> {
    get_expression_for_family_with_options(family, arguments, ExpressionOptions::default())
}

/// Enum-based alternative with explicit expansion resource limits.
pub fn get_expression_for_family_with_options(
    family: ScalarIntegral,
    arguments: &[Atom],
    options: ExpressionOptions,
) -> Result<LaurentSeries, String> {
    if arguments.len() != family.arity() {
        return Err(format!(
            "{} expects {} arguments including mu_squared last; received {}",
            family.name(),
            family.arity(),
            arguments.len()
        ));
    }
    if options.max_nodes == 0 || options.max_depth == 0 {
        return Err("expression max_nodes and max_depth must be positive".into());
    }
    let master = match family {
        ScalarIntegral::A0 => crate::A0(),
        ScalarIntegral::B0 => crate::B0(),
        ScalarIntegral::DB0 => crate::dB0(),
        ScalarIntegral::C0 => crate::C0(),
        ScalarIntegral::D0 => crate::D0(),
    };
    let mut expansion = Expansion::new(crate::expressions::shared_definitions(), options);
    let specialized = if family == ScalarIntegral::C0 {
        one_mass_one_scale_triangle(arguments)
    } else {
        None
    };
    let mut coefficients = Vec::with_capacity(3);
    for (index, tag) in [0, -1, -2].into_iter().enumerate() {
        let mut call = vec![Atom::num(tag)];
        call.extend_from_slice(arguments);
        let body = specialized.as_ref().map_or_else(
            || master.call(&call),
            |series| series.coefficients()[index].clone(),
        );
        let coefficient = expansion.expand(body.as_view(), 0, 0)?.atom;
        let mut leftover = None;
        coefficient.visitor(&mut |view| {
            if let Some(symbol) = view.get_symbol()
                && (symbol.get_name().contains("::__olo_")
                    || matches!(
                        symbol.get_name(),
                        "oneloopmaster::A0"
                            | "oneloopmaster::B0"
                            | "oneloopmaster::dB0"
                            | "oneloopmaster::C0"
                            | "oneloopmaster::D0"
                    ))
            {
                leftover = Some(symbol);
            }
            leftover.is_none()
        });
        if let Some(symbol) = leftover {
            return Err(format!(
                "incomplete native definition left unresolved symbol {symbol}"
            ));
        }
        coefficients.push(coefficient);
    }
    Ok(LaurentSeries::new(
        coefficients.remove(0),
        coefficients.remove(0),
        coefficients.remove(0),
    ))
}

/// Compact, integrated Feynman-parameter identity for the one-mass,
/// one-scale sector. This is an inspection-only specialization: the canonical
/// mapped definitions and all numerical backends are deliberately unchanged.
fn one_mass_one_scale_triangle(arguments: &[Atom]) -> Option<LaurentSeries> {
    let zero = Atom::num(0);
    let mass_index = (0..3).find(|&index| arguments[index + 3] != zero)?;
    if (0..3).any(|index| index != mass_index && arguments[index + 3] != zero) {
        return None;
    }
    let momentum_index = (mass_index + 1) % 3;
    if (0..3).any(|index| index != momentum_index && arguments[index] != zero) {
        return None;
    }
    let mass = &arguments[mass_index + 3];
    let momentum = &arguments[momentum_index];
    // No mass attribute is required: untyped masses may be genuinely complex.
    // As for all numeric entry points, the physical domain is Im(m²)<=0.
    if !known_real(momentum.as_view()) {
        return None;
    }
    let scale = &arguments[6];
    let pi_squared = Symbol::PI.to_atom().pow(2);
    // For Im(m)<0 and real s!=0, F=s*C0 obeys
    // F'=log(-s/m)/(m+s), F(s->0)=0, giving F=pi²/6-Li2(1+s/m).
    // Im(1+s/(m-i0)) has sign(s). Off the mass axis use ordinary
    // complex Li2; on it, select the upper lip for s>0 and lower for s<0.
    // The latter matters for negative real masses: always conjugating the
    // real-axis dilog would select the wrong sheet when both m and s are <0.
    let dilog = (Atom::num(1) + momentum / mass).polylog(2);
    let continued = crate::if_nonzero_else(
        &crate::sheet_exact::im(mass),
        dilog.clone(),
        crate::if_nonzero_else(
            &crate::sheet_exact::negative(momentum),
            dilog.clone(),
            dilog.conj(),
        ),
    );
    let finite = (pi_squared / 6 - continued) / momentum;
    let vacuum = LaurentSeries::new(
        (1 - crate::physical_log(&(mass / scale))) / mass,
        Atom::num(1) / mass,
        zero.clone(),
    );
    let massive =
        crate::triangle::if_series(momentum, crate::triangle::finite_series(finite), vacuum);
    let log = crate::physical_log(&(-momentum / scale));
    // These are already in the public Laurent normalization, including the
    // cancellation of the raw massless triangle's pi^2 double-pole correction.
    let massless = crate::triangle::if_series(
        momentum,
        LaurentSeries::new(
            log.pow(2) / (2 * momentum),
            -log / momentum,
            Atom::num(1) / momentum,
        ),
        crate::triangle::finite_series(zero),
    );
    Some(crate::triangle::if_series(mass, massive, massless))
}

impl ScalarIntegral {
    /// See [`get_expression`].
    pub fn get_expression(self, arguments: &[Atom]) -> Result<LaurentSeries, String> {
        get_expression_for_family(self, arguments)
    }
}

#[derive(Clone, Debug)]
struct Expanded {
    atom: Atom,
    nodes: usize,
}

fn helper_name(symbol: Symbol) -> Option<&'static str> {
    match symbol.get_name() {
        concat!(env!("CARGO_CRATE_NAME"), "::__olo_real_part") => Some("__olo_real_part"),
        concat!(env!("CARGO_CRATE_NAME"), "::__olo_imaginary_part") => Some("__olo_imaginary_part"),
        concat!(env!("CARGO_CRATE_NAME"), "::__olo_sign_nonnegative") => {
            Some("__olo_sign_nonnegative")
        }
        _ => None,
    }
}

fn known_real(value: AtomView<'_>) -> bool {
    value.is_real()
}

fn normalized_call(symbol: Symbol, mut arguments: Vec<Atom>) -> Atom {
    if arguments.len() == 1 {
        let argument = &arguments[0];
        if known_real(argument.as_view()) {
            if symbol == Symbol::CONJ || helper_name(symbol) == Some("__olo_real_part") {
                return arguments.remove(0);
            }
            if helper_name(symbol) == Some("__olo_imaginary_part") {
                return Atom::num(0);
            }
        }
        // Nonnegative is sufficient here: sign_nonnegative(0) is exactly +1.
        // This does not use positivity as a proof that a zero-test is true.
        if helper_name(symbol) == Some("__olo_sign_nonnegative") && argument.is_positive() {
            return Atom::num(1);
        }
    }
    symbol.call(&arguments)
}

// Arguments retain their caller's scope and are expanded only if referenced.
// This avoids expanding dead branches and preserves simultaneous substitution
// even when an actual argument happens to have a formal parameter's name.
struct Scope {
    parent: Option<usize>,
    bindings: Vec<(Atom, Atom, usize)>,
    definitions: Vec<Symbol>,
}

struct Expansion<'a> {
    map: &'a FunctionMap,
    options: ExpressionOptions,
    work: usize,
    scopes: Vec<Scope>,
    memo: HashMap<(usize, Atom), Expanded>,
    compact_memo: HashMap<(usize, Atom), Expanded>,
    calls: HashMap<Atom, Expanded>,
    closed_definitions: HashMap<Atom, bool>,
}

impl<'a> Expansion<'a> {
    fn new(map: &'a FunctionMap, options: ExpressionOptions) -> Self {
        Self {
            map,
            options,
            work: 0,
            scopes: vec![Scope {
                parent: None,
                bindings: vec![],
                definitions: vec![],
            }],
            memo: HashMap::new(),
            compact_memo: HashMap::new(),
            calls: HashMap::new(),
            closed_definitions: HashMap::new(),
        }
    }

    fn charge(&mut self, nodes: usize) -> Result<(), String> {
        self.work = self
            .work
            .checked_add(nodes)
            .ok_or_else(|| self.budget_error())?;
        if self.work > self.options.max_nodes {
            return Err(self.budget_error());
        }
        Ok(())
    }

    fn budget_error(&self) -> String {
        format!(
            "complete expression expansion exceeded max_nodes={}; raise the explicit budget or use compact OneLoopExpressions with its FunctionMap",
            self.options.max_nodes
        )
    }

    fn finish(&mut self, atom: Atom) -> Result<Expanded, String> {
        let mut nodes = 0usize;
        let remaining = self.options.max_nodes.saturating_sub(self.work);
        atom.visitor(&mut |_| {
            nodes = nodes.saturating_add(1);
            nodes <= remaining
        });
        self.charge(nodes)?;
        Ok(Expanded { atom, nodes })
    }

    fn binding(&self, atom: AtomView<'_>, mut scope: usize) -> Option<(Atom, usize)> {
        loop {
            let current = &self.scopes[scope];
            if let Some((_, value, caller)) =
                current.bindings.iter().find(|(key, _, _)| key == atom)
            {
                return Some((value.clone(), *caller));
            }
            scope = current.parent?;
        }
    }

    // A canonical explicit call is a complete memo key only when the body and
    // its callees do not capture parameters from an enclosing definition.
    fn closed_definition(&mut self, call: AtomView<'_>, active: &mut HashSet<Atom>) -> bool {
        let Some((tag_count, parameters, body)) = self.map.get_definition(call) else {
            return true;
        };
        let key = match call {
            AtomView::Fun(f) => f.get_symbol().call(
                &f.iter()
                    .take(tag_count)
                    .map(|a| a.to_owned())
                    .collect::<Vec<_>>(),
            ),
            _ => call.to_owned(),
        };
        if let Some(closed) = self.closed_definitions.get(&key) {
            return *closed;
        }
        if !active.insert(key.clone()) {
            return false;
        }
        let mut closed = true;
        let mut callees = Vec::new();
        body.visitor(&mut |view| {
            if parameters
                .iter()
                .any(|parameter| parameter.as_view() == view)
            {
                return false;
            }
            if let AtomView::Var(variable) = view
                && !variable.get_symbol().is_builtin()
                && self.map.get_definition(view).is_none()
            {
                closed = false;
            }
            if self.map.get_definition(view).is_some() {
                callees.push(view.to_owned());
            }
            closed
        });
        if closed {
            for callee in callees {
                if !self.closed_definition(callee.as_view(), active) {
                    closed = false;
                    break;
                }
            }
        }
        active.remove(&key);
        self.closed_definitions.insert(key, closed);
        closed
    }

    // Substitute actual arguments without expanding registered definitions.
    // Equivalent calls from different lexical scopes then share one result.
    // Numeric IF conditions still prune their dead arms before substitution.
    fn compact(
        &mut self,
        input: AtomView<'_>,
        scope: usize,
        depth: usize,
    ) -> Result<Expanded, String> {
        if depth > self.options.max_depth {
            return Err(format!(
                "complete expression expansion exceeded max_depth={}",
                self.options.max_depth
            ));
        }
        self.charge(1)?;
        let key = (scope, input.to_owned());
        if let Some(cached) = self.compact_memo.get(&key) {
            let cached = cached.clone();
            self.charge(cached.nodes)?;
            return Ok(cached);
        }
        let result = if let Some((value, caller)) = self.binding(input, scope) {
            self.compact(value.as_view(), caller, depth + 1)?
        } else {
            match input {
                AtomView::Num(_) | AtomView::Var(_) => self.finish(input.to_owned())?,
                AtomView::Fun(f) if f.get_symbol() == Symbol::IF && f.get_nargs() == 3 => {
                    let args = f.iter().collect::<Vec<_>>();
                    let condition = self.compact(args[0], scope, depth + 1)?;
                    if let AtomView::Num(number) = condition.atom.as_view() {
                        self.compact(args[if number.is_zero() { 2 } else { 1 }], scope, depth + 1)?
                    } else {
                        let yes = self.compact(args[1], scope, depth + 1)?;
                        let no = self.compact(args[2], scope, depth + 1)?;
                        self.finish(Symbol::IF.call((condition.atom, yes.atom, no.atom)))?
                    }
                }
                AtomView::Fun(f) => {
                    let mut arguments = Vec::with_capacity(f.get_nargs());
                    for argument in f.iter() {
                        arguments.push(self.compact(argument, scope, depth + 1)?.atom);
                    }
                    let result = normalized_call(f.get_symbol(), arguments);
                    self.finish(result)?
                }
                AtomView::Pow(p) => {
                    let (base, exponent) = p.get_base_exp();
                    let base = self.compact(base, scope, depth + 1)?;
                    let exponent = self.compact(exponent, scope, depth + 1)?;
                    self.finish(base.atom.pow(exponent.atom))?
                }
                AtomView::Mul(m) => {
                    let mut values = Vec::new();
                    for child in m.iter() {
                        values.push(self.compact(child, scope, depth + 1)?.atom);
                    }
                    self.finish(values.into_iter().product())?
                }
                AtomView::Add(a) => {
                    let mut values = Vec::new();
                    for child in a.iter() {
                        values.push(self.compact(child, scope, depth + 1)?.atom);
                    }
                    self.finish(values.into_iter().sum())?
                }
            }
        };
        self.compact_memo.insert(key, result.clone());
        Ok(result)
    }

    fn expand(
        &mut self,
        input: AtomView<'_>,
        scope: usize,
        depth: usize,
    ) -> Result<Expanded, String> {
        if depth > self.options.max_depth {
            return Err(format!(
                "complete expression expansion exceeded max_depth={}",
                self.options.max_depth
            ));
        }
        self.charge(1)?;
        let key = (scope, input.to_owned());
        if let Some(cached) = self.memo.get(&key) {
            let cached = cached.clone();
            self.charge(cached.nodes)?;
            return Ok(cached);
        }
        let result = if let Some((value, caller)) = self.binding(input, scope) {
            self.expand(value.as_view(), caller, depth + 1)?
        } else if let Some((tag_count, parameters, body)) = self.map.get_definition(input) {
            let symbol = input.get_symbol().expect("a definition has a head symbol");
            if self.scopes[scope].definitions.contains(&symbol) {
                return Err(format!("cyclic native definition encountered at {symbol}"));
            }
            let actual = match input {
                AtomView::Fun(f) => f
                    .iter()
                    .skip(tag_count)
                    .map(|a| a.to_owned())
                    .collect::<Vec<_>>(),
                _ => vec![],
            };
            if actual.len() != parameters.len() {
                return Err(format!(
                    "native definition {symbol} expects {} untagged arguments; received {}",
                    parameters.len(),
                    actual.len()
                ));
            }
            let body = body.clone();
            let mut bindings: Vec<_> = parameters
                .iter()
                .zip(actual)
                .map(|(parameter, value)| (parameter.as_view().to_owned(), value, scope))
                .collect();
            let closed = self.closed_definition(input, &mut HashSet::new());
            let canonical = if closed {
                Some(self.compact(input, scope, depth + 1)?.atom)
            } else {
                None
            };
            if let Some(cached) = canonical.as_ref().and_then(|key| self.calls.get(key)) {
                let cached = cached.clone();
                self.charge(cached.nodes)?;
                self.memo.insert(key, cached.clone());
                return Ok(cached);
            }
            if let Some(canonical) = &canonical {
                if canonical.as_view().get_symbol() != Some(symbol) {
                    let result = self.expand(canonical.as_view(), 0, depth + 1)?;
                    self.calls.insert(canonical.clone(), result.clone());
                    self.memo.insert(key, result.clone());
                    return Ok(result);
                }
                if let AtomView::Fun(call) = canonical.as_view() {
                    for ((_, actual, caller), canonical_argument) in
                        bindings.iter_mut().zip(call.iter().skip(tag_count))
                    {
                        *actual = canonical_argument.to_owned();
                        *caller = 0;
                    }
                }
            }
            let mut definitions = self.scopes[scope].definitions.clone();
            definitions.push(symbol);
            let child = self.scopes.len();
            self.scopes.push(Scope {
                parent: Some(scope),
                bindings,
                definitions,
            });
            let result = self.expand(body.as_view(), child, depth + 1)?;
            if let Some(canonical) = canonical {
                self.calls.insert(canonical, result.clone());
            }
            result
        } else {
            match input {
                AtomView::Num(_) | AtomView::Var(_) => self.finish(input.to_owned())?,
                AtomView::Fun(f) if f.get_symbol() == Symbol::IF && f.get_nargs() == 3 => {
                    let arguments = f.iter().collect::<Vec<_>>();
                    let condition = self.expand(arguments[0], scope, depth + 1)?;
                    if let AtomView::Num(number) = condition.atom.as_view() {
                        // Native IF tests exact zero. Positive attributes are NOT
                        // a nonzero proof: Symbolica's positivity includes zero.
                        let chosen = if number.is_zero() { 2 } else { 1 };
                        self.expand(arguments[chosen], scope, depth + 1)?
                    } else {
                        let yes = self.expand(arguments[1], scope, depth + 1)?;
                        let no = self.expand(arguments[2], scope, depth + 1)?;
                        self.finish(Symbol::IF.call((condition.atom, yes.atom, no.atom)))?
                    }
                }
                AtomView::Fun(f) => {
                    if f.get_symbol().get_name().contains("::__olo_") {
                        return Err(format!(
                            "no native definition for helper {}",
                            f.get_symbol()
                        ));
                    }
                    let mut arguments = Vec::with_capacity(f.get_nargs());
                    for argument in f.iter() {
                        arguments.push(self.expand(argument, scope, depth + 1)?.atom);
                    }
                    let result = normalized_call(f.get_symbol(), arguments);
                    self.finish(result)?
                }
                AtomView::Pow(p) => {
                    let (base, exponent) = p.get_base_exp();
                    let base = self.expand(base, scope, depth + 1)?;
                    let exponent = self.expand(exponent, scope, depth + 1)?;
                    self.finish(base.atom.pow(exponent.atom))?
                }
                AtomView::Mul(m) => {
                    let mut values = Vec::new();
                    for child in m.iter() {
                        values.push(self.expand(child, scope, depth + 1)?.atom);
                    }
                    self.finish(values.into_iter().product())?
                }
                AtomView::Add(a) => {
                    let mut values = Vec::new();
                    for child in a.iter() {
                        values.push(self.expand(child, scope, depth + 1)?.atom);
                    }
                    self.finish(values.into_iter().sum())?
                }
            }
        };
        self.memo.insert(key, result.clone());
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substitution_is_scoped_and_dead_branches_are_not_expanded() {
        let x = symbolica::symbol!("inspection_scope_x");
        let y = symbolica::symbol!("inspection_scope_y");
        let f = symbolica::symbol!("inspection_scope_f");
        let cycle = symbolica::symbol!("inspection_scope_cycle");
        let mut map = FunctionMap::new();
        map.add_function(f, vec![x, y], x.to_atom() - y).unwrap();
        map.add_function(cycle, vec![x], cycle.call(x)).unwrap();
        let mut expansion = Expansion::new(&map, ExpressionOptions::default());
        assert_eq!(
            expansion
                .expand(f.call((y, x)).as_view(), 0, 0)
                .unwrap()
                .atom,
            y.to_atom() - x
        );
        let nested = f.call((f.call((3, 1)), 1));
        assert_eq!(expansion.expand(nested.as_view(), 0, 0).unwrap().atom, 1);
        let branch = Symbol::IF.call((x, cycle.call(x), 7));
        let g = symbolica::symbol!("inspection_scope_g");
        map.add_function(g, vec![x], branch).unwrap();
        let mut expansion = Expansion::new(&map, ExpressionOptions::default());
        assert_eq!(expansion.expand(g.call(0).as_view(), 0, 0).unwrap().atom, 7);
        assert!(
            expansion
                .expand(g.call(1).as_view(), 0, 0)
                .unwrap_err()
                .contains("cyclic")
        );
        // A FunctionMap helper may capture its caller's formal parameters.
        // Canonical-call memoization must not confuse those distinct scopes.
        let captured = symbolica::symbol!("inspection_scope_captured");
        let indirect = symbolica::symbol!("inspection_scope_indirect");
        let caller = symbolica::symbol!("inspection_scope_caller");
        map.add_function(captured, Vec::<Symbol>::new(), x.to_atom())
            .unwrap();
        map.add_function(
            indirect,
            Vec::<Symbol>::new(),
            captured.call(&[] as &[Atom]),
        )
        .unwrap();
        map.add_function(caller, vec![x], indirect.call(&[] as &[Atom]))
            .unwrap();
        let mut expansion = Expansion::new(&map, ExpressionOptions::default());
        for value in [1, 2] {
            assert_eq!(
                expansion
                    .expand(caller.call(value).as_view(), 0, 0)
                    .unwrap()
                    .atom,
                value
            );
        }
    }
}
