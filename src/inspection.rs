//! Complete, inspectable master expressions, independent of numerical backends.
use crate::definitions::FunctionMap;
use crate::{LaurentSeries, ScalarIntegral};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    id::Replacement,
    transcendental::TranscendentalFunctions,
};

#[path = "inspection_graph.rs"]
mod graph;
pub use graph::{
    SharedLaurentSeries, get_expression_shared, get_expression_shared_with_options,
    select_branch_shared,
};

pub(crate) fn prepare_symbols() {
    let _ = crate::inspection_formulas::definitions();
    let _ = graph::alias_symbols();
    let _ = symbolica::symbol!("oneloopmaster::__olo_expression_thunk");
}

/// Resource limits for fully substituting the transparent native definitions.
///
/// Generic boxes can be much larger than their compact FunctionMap form. Limits
/// fail explicitly; they never return a truncated formula or opaque helper.
#[derive(Clone, Copy, Debug)]
pub struct ExpressionOptions {
    /// Maximum nodes in a materialized expression and unique expansion work.
    /// Reusing a cached subtree does not consume its entire size again.
    pub max_nodes: usize,
    /// Maximum recursive expression/definition depth.
    pub max_depth: usize,
}

impl Default for ExpressionOptions {
    fn default() -> Self {
        Self {
            max_nodes: 10_000_000,
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

/// Expand a master only on the region selected by condition-only probe rules.
///
/// Unlike selecting after [`get_expression`], discarded branches are never
/// inlined. The rules have the same semantics as [`crate::select_branch`]: the
/// returned bodies remain parametric and undecidable conditions remain `if`s.
pub fn get_expression_on_branch(
    master: impl AtomCore,
    branch_rules: &[Replacement],
) -> Result<LaurentSeries, String> {
    get_expression_on_branch_with_options(master, branch_rules, ExpressionOptions::default())
}

/// Version of [`get_expression_on_branch`] with explicit expansion limits.
pub fn get_expression_on_branch_with_options(
    master: impl AtomCore,
    branch_rules: &[Replacement],
    options: ExpressionOptions,
) -> Result<LaurentSeries, String> {
    let (family, arguments) = master_arguments(master)?;
    expand_family(family, &arguments, options, Some(branch_rules))
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
    expand_family(family, arguments, options, None)
}

fn expand_family(
    family: ScalarIntegral,
    arguments: &[Atom],
    options: ExpressionOptions,
    branch_rules: Option<&[Replacement]>,
) -> Result<LaurentSeries, String> {
    expand_family_impl(family, arguments, options, branch_rules, None).map(|(series, _)| series)
}

fn expand_family_impl(
    family: ScalarIntegral,
    arguments: &[Atom],
    options: ExpressionOptions,
    branch_rules: Option<&[Replacement]>,
    graph: Option<graph::GraphBuilder>,
) -> Result<(LaurentSeries, Option<graph::GraphBuilder>), String> {
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
    let mut expansion = Expansion::new(crate::inspection_formulas::definitions(), options);
    expansion.graph = graph;
    expansion.branch_rules = branch_rules.map(|rules| {
        rules
            .iter()
            .map(crate::branch_selection::exact_replacement)
            .collect()
    });
    let specialized = if family == ScalarIntegral::C0 {
        opposite_invariants_triangle(arguments).or_else(|| one_mass_one_scale_triangle(arguments))
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
        let coefficient = expansion
            .expand(body.as_view(), 0, 0)?
            .atom
            .as_ref()
            .clone();
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
    Ok((
        LaurentSeries::new(
            coefficients.remove(0),
            coefficients.remove(0),
            coefficients.remove(0),
        ),
        expansion.graph,
    ))
}

/// Integrate the linear Feynman parameter first for C0(0,-a,a,a,a,b).
/// The remaining logarithm factors into two quadratics, giving four dilogs.
/// This includes the a=0 sector and both real-axis lips; it is not a probe-
/// selected formula. For physical masses Im(a),Im(b)<=0, b/a lies off the
/// negative real cut unless both masses are real and have opposite signs.
fn opposite_invariants_triangle(arguments: &[Atom]) -> Option<LaurentSeries> {
    let a = &arguments[2];
    let b = &arguments[5];
    if !arguments[0].is_zero() || arguments[1] != -a || arguments[3] != a || arguments[4] != a {
        return None;
    }
    let r = b / a;
    let d = (r.pow(2) + 4).sqrt();
    let e = (&r * (&r - 4)).sqrt();
    // On r<0 each dilog cut argument approaches from the side opposite to r.
    // For real physical masses that side is upper for a>0, lower for a<0.
    let lip = crate::sheet_exact::sign_nonnegative(&crate::sheet_exact::re(a));
    let li2 = |z: Atom| {
        let cut =
            Symbol::PI.to_atom().pow(2) / 6 - crate::dilog_atom(1 - &z) - z.log() * (&z - 1).log()
                + crate::sheet_exact::i() * Symbol::PI.to_atom() * &lip * z.log();
        let real = crate::if_nonzero_else(
            &crate::sheet_exact::negative(&(1 - &z)),
            cut,
            crate::dilog_atom(&z),
        );
        crate::if_nonzero_else(&crate::sheet_exact::im(&z), crate::dilog_atom(&z), real)
    };
    let finite =
        (li2((-&r + &d) / 2) + li2((-&r - d) / 2) - li2((2 - &r + &e) / 2) - li2((2 - &r - e) / 2))
            / (2 * a);
    let zero_a_finite = crate::if_nonzero(b, (1 - crate::physical_log(&(b / &arguments[6]))) / b);
    Some(LaurentSeries::new(
        crate::if_nonzero_else(a, finite, zero_a_finite),
        crate::if_nonzero_else(a, Atom::num(0), crate::if_nonzero(b, Atom::num(1) / b)),
        Atom::num(0),
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
    atom: Arc<Atom>,
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
    value.is_real().is_true()
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
        if helper_name(symbol) == Some("__olo_sign_nonnegative") && argument.is_positive().is_true()
        {
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
    memo: HashMap<(usize, usize, bool, Atom), Expanded>,
    compact_memo: HashMap<(usize, usize, bool, Atom), Expanded>,
    calls: HashMap<(usize, bool, Arc<Atom>), Expanded>,
    closed_definitions: HashMap<Atom, bool>,
    used_parameters: HashMap<Atom, Vec<bool>>,
    // A fresh expansion owns its probe context. No cache is shared across calls
    // with different rules or between all-branches and branch-selected modes.
    branch_rules: Option<Vec<Replacement>>,
    condition_memo: HashMap<Atom, Option<bool>>,
    graph: Option<graph::GraphBuilder>,
    aliases: HashMap<Atom, Atom>,
    active_aliases: HashSet<Atom>,
    context: usize,
    contexts: Vec<Option<(usize, Atom, bool)>>,
    compact_thunks: HashMap<Atom, Atom>,
    thunk_intern: HashMap<Atom, Atom>,
    preserve_predicate: bool,
    graph_probe: Option<graph::NumericProbe>,
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
            used_parameters: HashMap::new(),
            branch_rules: None,
            condition_memo: HashMap::new(),
            graph: None,
            aliases: HashMap::new(),
            active_aliases: HashSet::new(),
            context: 0,
            contexts: vec![None],
            compact_thunks: HashMap::new(),
            thunk_intern: HashMap::new(),
            preserve_predicate: false,
            graph_probe: None,
        }
    }

    fn condition_truth(&mut self, condition: &Atom) -> Option<bool> {
        if self.preserve_predicate {
            return if let AtomView::Num(number) = condition.as_view() {
                Some(!number.is_zero())
            } else {
                None
            };
        }
        let mut context = self.context;
        while let Some((parent, predicate, truth)) = &self.contexts[context] {
            if predicate == condition {
                return Some(*truth);
            }
            context = *parent;
        }
        if let Some(rules) = &self.branch_rules {
            if let Some(value) = self.condition_memo.get(condition) {
                return *value;
            }
            let sample = condition.replace_multiple(rules);
            let value = crate::branch_selection::numeric_truth(sample.as_view());
            self.condition_memo.insert(condition.clone(), value);
            value
        } else if let AtomView::Num(number) = condition.as_view() {
            Some(!number.is_zero())
        } else {
            None
        }
    }

    fn expand_arm(
        &mut self,
        input: AtomView<'_>,
        scope: usize,
        depth: usize,
        condition: &Atom,
        truth: bool,
    ) -> Result<Expanded, String> {
        if self.graph.is_some() || self.preserve_predicate {
            return self.expand(input, scope, depth);
        }
        let previous = self.context;
        self.context = self.contexts.len();
        self.contexts
            .push(Some((previous, condition.clone(), truth)));
        let result = self.expand(input, scope, depth);
        self.context = previous;
        result
    }

    fn expand_condition(
        &mut self,
        input: AtomView<'_>,
        scope: usize,
        depth: usize,
    ) -> Result<Expanded, String> {
        let previous = self.preserve_predicate;
        self.preserve_predicate |= self.branch_rules.is_some();
        let result = self.expand(input, scope, depth);
        self.preserve_predicate = previous;
        result
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
            "complete expression expansion exceeded max_nodes={} ({} expansion visits); raise the explicit budget or select a branch during construction",
            self.options.max_nodes, self.work
        )
    }

    fn finish(&mut self, atom: Atom) -> Result<Expanded, String> {
        // Serialized bytes bound the node count from above. Only traverse a
        // large result when the inexpensive byte bound cannot certify it.
        // In particular, do not recount every child whenever a parent is built.
        if atom.as_view().get_byte_size() > self.options.max_nodes {
            let mut nodes = 0usize;
            atom.visitor(&mut |_| {
                nodes = nodes.saturating_add(1);
                nodes <= self.options.max_nodes
            });
            if nodes > self.options.max_nodes {
                return Err(format!(
                    "{}; materialized subtree has at least {nodes} nodes ({} bytes)",
                    self.budget_error(),
                    atom.as_view().get_byte_size()
                ));
            }
        }
        Ok(Expanded {
            atom: Arc::new(atom),
        })
    }

    fn finish_expanded(&mut self, atom: Atom) -> Result<Expanded, String> {
        let atom = if let Some(graph) = &mut self.graph {
            graph.intern(atom)
        } else {
            atom
        };
        self.finish(atom)
    }

    fn finish_compact(&mut self, atom: Atom) -> Result<Expanded, String> {
        if let Some(graph) = &self.graph
            && !matches!(atom.as_view(), AtomView::Num(_) | AtomView::Var(_))
            && !self.compact_thunks.contains_key(&atom)
        {
            if let Some(alias) = self.thunk_intern.get(&atom) {
                return Ok(Expanded {
                    atom: Arc::new(alias.clone()),
                });
            }
            let alias = symbolica::symbol!("oneloopmaster::__olo_expression_thunk")
                .call((graph.id, self.compact_thunks.len()));
            self.compact_thunks.insert(alias.clone(), atom.clone());
            self.thunk_intern.insert(atom, alias.clone());
            return Ok(Expanded {
                atom: Arc::new(alias),
            });
        }
        self.finish(atom)
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

    // A sheet component often uses only two of its six formal arguments.
    // Do not copy the other four (recursively growing sheet histories) merely
    // to construct a memo key. This changes neither the body nor its branches.
    fn used_parameters(
        &mut self,
        call: AtomView<'_>,
        active: &mut HashSet<Atom>,
    ) -> Option<Vec<bool>> {
        // Open definitions may capture a caller's parameter through a callee
        // with no explicit arguments. Keep all arguments in that case.
        if !self.closed_definition(call, &mut HashSet::new()) {
            return None;
        }
        let (tags, parameters, body) = self.map.get_definition(call)?;
        let key = match call {
            AtomView::Fun(f) => f.get_symbol().call(
                &f.iter()
                    .take(tags)
                    .map(|a| a.to_owned())
                    .collect::<Vec<_>>(),
            ),
            _ => call.to_owned(),
        };
        if let Some(used) = self.used_parameters.get(&key) {
            return Some(used.clone());
        }
        if !active.insert(key.clone()) {
            return Some(vec![true; parameters.len()]);
        }
        let parameters = parameters
            .iter()
            .map(|p| p.as_view().to_owned())
            .collect::<Vec<_>>();
        let body = body.clone();
        let mut used = vec![false; parameters.len()];
        let mut pending = vec![body.as_view()];
        while let Some(view) = pending.pop() {
            if let Some(i) = parameters.iter().position(|p| p.as_view() == view) {
                used[i] = true;
                continue;
            }
            match view {
                AtomView::Fun(f) => {
                    if let Some(callee_used) = self.used_parameters(view, active) {
                        let (tags, _, _) = self.map.get_definition(view).unwrap();
                        pending.extend(f.iter().take(tags));
                        pending.extend(
                            f.iter()
                                .skip(tags)
                                .zip(callee_used)
                                .filter_map(|(a, used)| used.then_some(a)),
                        );
                    } else {
                        pending.extend(f.iter());
                    }
                }
                AtomView::Pow(p) => {
                    let (a, b) = p.get_base_exp();
                    pending.extend([a, b]);
                }
                AtomView::Mul(m) => pending.extend(m.iter()),
                AtomView::Add(a) => pending.extend(a.iter()),
                _ => {}
            }
        }
        active.remove(&key);
        self.used_parameters.insert(key, used.clone());
        Some(used)
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
        let key = (
            scope,
            self.context,
            self.preserve_predicate,
            input.to_owned(),
        );
        if let Some(cached) = self.compact_memo.get(&key) {
            let cached = cached.clone();
            return Ok(cached);
        }
        let result = if self.compact_thunks.contains_key(&input.to_owned()) {
            self.finish(input.to_owned())?
        } else if let Some((value, caller)) = self.binding(input, scope) {
            self.compact(value.as_view(), caller, depth + 1)?
        } else {
            match input {
                AtomView::Num(_) | AtomView::Var(_) => self.finish_compact(input.to_owned())?,
                AtomView::Fun(f) if f.get_symbol() == Symbol::IF && f.get_nargs() == 3 => {
                    let args = f.iter().collect::<Vec<_>>();
                    let condition = self.compact(args[0], scope, depth + 1)?;
                    if let Some(nonzero) = self.condition_truth(&condition.atom) {
                        self.compact(args[if nonzero { 1 } else { 2 }], scope, depth + 1)?
                    } else {
                        let yes = self.compact(args[1], scope, depth + 1)?;
                        let no = self.compact(args[2], scope, depth + 1)?;
                        self.finish_compact(Symbol::IF.call((
                            condition.atom.as_ref(),
                            yes.atom.as_ref(),
                            no.atom.as_ref(),
                        )))?
                    }
                }
                AtomView::Fun(f) => {
                    let mut arguments = Vec::with_capacity(f.get_nargs());
                    let used = self.used_parameters(input, &mut HashSet::new());
                    let tags = self.map.get_definition(input).map_or(0, |d| d.0);
                    for (i, argument) in f.iter().enumerate() {
                        if i >= tags && used.as_ref().is_some_and(|u| !u[i - tags]) {
                            arguments.push(Arc::new(Atom::num(0)));
                        } else {
                            arguments.push(self.compact(argument, scope, depth + 1)?.atom);
                        }
                    }
                    let result = normalized_call(
                        f.get_symbol(),
                        arguments.iter().map(|a| a.as_ref().clone()).collect(),
                    );
                    self.finish_compact(result)?
                }
                AtomView::Pow(p) => {
                    let (base, exponent) = p.get_base_exp();
                    let base = self.compact(base, scope, depth + 1)?;
                    let exponent = self.compact(exponent, scope, depth + 1)?;
                    self.finish_compact(base.atom.pow(exponent.atom.as_ref()))?
                }
                AtomView::Mul(m) => {
                    let mut values = Vec::new();
                    for child in m.iter() {
                        values.push(self.compact(child, scope, depth + 1)?.atom);
                    }
                    self.finish_compact(values.iter().map(|a| a.as_ref()).product())?
                }
                AtomView::Add(a) => {
                    let mut values = Vec::new();
                    for child in a.iter() {
                        values.push(self.compact(child, scope, depth + 1)?.atom);
                    }
                    self.finish_compact(values.iter().map(|a| a.as_ref()).sum())?
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
        let key = (
            scope,
            self.context,
            self.preserve_predicate,
            input.to_owned(),
        );
        if let Some(cached) = self.memo.get(&key) {
            let cached = cached.clone();
            return Ok(cached);
        }
        let result = if let Some((value, caller)) = self.binding(input, scope) {
            self.expand(value.as_view(), caller, depth + 1)?
        } else if let Some(body) = self
            .aliases
            .get(&input.to_owned())
            .or_else(|| self.compact_thunks.get(&input.to_owned()))
            .cloned()
        {
            let alias = input.to_owned();
            if !self.active_aliases.insert(alias.clone()) {
                return Err(format!("cyclic expression alias {alias}"));
            }
            let result = self.expand(body.as_view(), 0, depth + 1)?;
            self.active_aliases.remove(&alias);
            result
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
                let mut canonical = self.compact(input, scope, depth + 1)?.atom;
                while let Some(body) = self.compact_thunks.get(canonical.as_ref()) {
                    canonical = Arc::new(body.clone());
                }
                Some(canonical)
            } else {
                None
            };
            if let Some(cached) = canonical.as_ref().and_then(|key| {
                self.calls
                    .get(&(self.context, self.preserve_predicate, key.clone()))
            }) {
                let cached = cached.clone();
                self.memo.insert(key, cached.clone());
                return Ok(cached);
            }
            if let Some(canonical) = &canonical {
                if canonical.as_view().get_symbol() != Some(symbol) {
                    let result = self.expand(canonical.as_view(), 0, depth + 1)?;
                    self.calls.insert(
                        (self.context, self.preserve_predicate, canonical.clone()),
                        result.clone(),
                    );
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
            let result = self
                .expand(body.as_view(), child, depth + 1)
                .map_err(|error| format!("{error}\n  while expanding {symbol}"))?;
            if let Some(canonical) = canonical {
                self.calls.insert(
                    (self.context, self.preserve_predicate, canonical),
                    result.clone(),
                );
            }
            result
        } else {
            match input {
                AtomView::Num(_) | AtomView::Var(_) => self.finish_expanded(input.to_owned())?,
                AtomView::Fun(f) if f.get_symbol() == Symbol::IF && f.get_nargs() == 3 => {
                    let arguments = f.iter().collect::<Vec<_>>();
                    if !self.preserve_predicate
                        && let Some(probe) = &mut self.graph_probe
                        && let Some(truth) = probe.truth(
                            arguments[0],
                            &self.aliases,
                            self.branch_rules.as_ref().unwrap(),
                        )
                    {
                        let result =
                            self.expand(arguments[if truth { 1 } else { 2 }], scope, depth + 1)?;
                        self.memo.insert(key, result.clone());
                        return Ok(result);
                    }
                    let condition = self.expand_condition(arguments[0], scope, depth + 1)?;
                    if let Some(nonzero) = self.condition_truth(&condition.atom) {
                        // Native IF tests exact zero. Positive attributes are NOT
                        // a nonzero proof: Symbolica's positivity includes zero.
                        let chosen = if nonzero { 1 } else { 2 };
                        self.expand(arguments[chosen], scope, depth + 1)?
                    } else {
                        let yes =
                            self.expand_arm(arguments[1], scope, depth + 1, &condition.atom, true)?;
                        let no = self.expand_arm(
                            arguments[2],
                            scope,
                            depth + 1,
                            &condition.atom,
                            false,
                        )?;
                        self.finish_expanded(Symbol::IF.call((
                            condition.atom.as_ref(),
                            yes.atom.as_ref(),
                            no.atom.as_ref(),
                        )))?
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
                    let result = normalized_call(
                        f.get_symbol(),
                        arguments.iter().map(|a| a.as_ref().clone()).collect(),
                    );
                    self.finish_expanded(result)?
                }
                AtomView::Pow(p) => {
                    let (base, exponent) = p.get_base_exp();
                    let base = self.expand(base, scope, depth + 1)?;
                    let exponent = self.expand(exponent, scope, depth + 1)?;
                    self.finish_expanded(base.atom.pow(exponent.atom.as_ref()))?
                }
                AtomView::Mul(m) => {
                    let mut values = Vec::new();
                    for child in m.iter() {
                        values.push(self.expand(child, scope, depth + 1)?.atom);
                    }
                    self.finish_expanded(values.iter().map(|a| a.as_ref()).product())?
                }
                AtomView::Add(a) => {
                    let mut values = Vec::new();
                    for child in a.iter() {
                        values.push(self.expand(child, scope, depth + 1)?.atom);
                    }
                    self.finish_expanded(values.iter().map(|a| a.as_ref()).sum())?
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
    fn probe_rules_prune_dead_definitions_and_never_substitute_bodies() {
        let x = symbolica::symbol!("inspection_probe_x");
        let y = symbolica::symbol!("inspection_probe_y");
        let f = symbolica::symbol!("inspection_probe_f");
        let cycle = symbolica::symbol!("inspection_probe_cycle");
        let mut map = FunctionMap::new();
        map.add_function(cycle, vec![x], cycle.call(x)).unwrap();
        let nested = Symbol::IF.call((y, x.to_atom() + y, x));
        map.add_function(f, vec![x], Symbol::IF.call((x, nested, cycle.call(x))))
            .unwrap();
        let rules = [
            Replacement::new(x.to_atom(), 2),
            Replacement::new(y.to_atom(), 3),
        ];
        let mut expansion = Expansion::new(
            &map,
            ExpressionOptions {
                max_nodes: 1_000,
                max_depth: 20,
            },
        );
        expansion.branch_rules = Some(
            rules
                .iter()
                .map(crate::branch_selection::exact_replacement)
                .collect(),
        );
        assert_eq!(
            expansion
                .expand(f.call(x).as_view(), 0, 0)
                .unwrap()
                .atom
                .as_ref()
                .clone(),
            x.to_atom() + y
        );
        // A different probe has a distinct cache and correctly exposes a cycle.
        let mut expansion = Expansion::new(&map, ExpressionOptions::default());
        expansion.branch_rules = Some(vec![crate::branch_selection::exact_replacement(
            &Replacement::new(x.to_atom(), 0),
        )]);
        assert!(
            expansion
                .expand(f.call(x).as_view(), 0, 0)
                .unwrap_err()
                .contains("cyclic")
        );
        // Unresolved conditions keep the original symbolic expression, while
        // nested decidable IFs still disappear.
        let expression = Symbol::IF.call((x.to_atom() + y, Symbol::IF.call((x, y, 17)), x));
        assert_eq!(
            expansion
                .expand(expression.as_view(), 0, 0)
                .unwrap()
                .atom
                .as_ref()
                .clone(),
            Symbol::IF.call((x.to_atom() + y, 17, x))
        );
    }

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
                .atom
                .as_ref()
                .clone(),
            y.to_atom() - x
        );
        let nested = f.call((f.call((3, 1)), 1));
        assert_eq!(
            expansion
                .expand(nested.as_view(), 0, 0)
                .unwrap()
                .atom
                .as_ref(),
            &Atom::num(1)
        );
        let branch = Symbol::IF.call((x, cycle.call(x), 7));
        let g = symbolica::symbol!("inspection_scope_g");
        map.add_function(g, vec![x], branch).unwrap();
        let mut expansion = Expansion::new(&map, ExpressionOptions::default());
        assert_eq!(
            expansion
                .expand(g.call(0).as_view(), 0, 0)
                .unwrap()
                .atom
                .as_ref(),
            &Atom::num(7)
        );
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
                    .atom
                    .as_ref()
                    .clone(),
                value
            );
        }
    }
}
