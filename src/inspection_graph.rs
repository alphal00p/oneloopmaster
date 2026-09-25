//! Complete expressions with native Symbolica common-subexpression bindings.
use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::atom::AliasedAtom;
use symbolica::domains::float::{Complex, Float, SingleFloat};

pub(super) fn alias_symbols() -> &'static [Symbol; 3] {
    use symbolica::atom::{SymbolAttribute, SymbolBuilder};
    static SYMBOLS: std::sync::OnceLock<[Symbol; 3]> = std::sync::OnceLock::new();
    SYMBOLS.get_or_init(|| {
        [
            ("oneloopmaster::t", vec![]),
            ("oneloopmaster::t_real", vec![SymbolAttribute::Real]),
            ("oneloopmaster::t_positive", vec![SymbolAttribute::Positive]),
        ]
        .map(|(name, attributes)| {
            SymbolBuilder::new(symbolica::wrap_symbol!(name))
                .with_attributes(attributes)
                .with_print_function(|view, options, _| {
                    if options.include_attributes {
                        return None;
                    }
                    let AtomView::Fun(f) = view else {
                        return None;
                    };
                    if f.get_nargs() != 2 {
                        return None;
                    }
                    let mut args = f.iter();
                    let id = args.next()?;
                    let index = args.next()?;
                    if options.mode.is_latex() {
                        Some(format!("t^{{({id})}}_{{{index}}}"))
                    } else if options.mode.is_symbolica() {
                        Some(format!("t[{id},{index}]"))
                    } else {
                        None
                    }
                })
                .build()
                .expect("unique transparent expression reference symbol")
        })
    })
}

/// All three complete coefficients, stored as Symbolica `AliasedAtom`s.
///
/// Every binding is included and contains only native expressions and earlier
/// bindings. These are not master calls, evaluation hooks, or opaque helpers.
/// Use the native `.evaluator(...)` on a coefficient, inspect `get_aliases()`,
/// or select its branches with [`select_branch_shared`] to obtain a plain Atom.
#[derive(Clone, Debug)]
pub struct SharedLaurentSeries {
    coefficients: [AliasedAtom; 3],
}

impl SharedLaurentSeries {
    pub fn coefficients(&self) -> &[AliasedAtom; 3] {
        &self.coefficients
    }
    pub fn into_coefficients(self) -> [AliasedAtom; 3] {
        self.coefficients
    }
}

/// Construct every branch without duplicating repeated subexpressions.
pub fn get_expression_shared(master: impl AtomCore) -> Result<SharedLaurentSeries, String> {
    get_expression_shared_with_options(master, ExpressionOptions::default())
}

/// Version of [`get_expression_shared`] with explicit resource limits.
/// The node budget concerns the shared construction, not a duplicated tree.
pub fn get_expression_shared_with_options(
    master: impl AtomCore,
    options: ExpressionOptions,
) -> Result<SharedLaurentSeries, String> {
    let (family, arguments) = master_arguments(master)?;
    let (series, graph) =
        expand_family_impl(family, &arguments, options, None, Some(GraphBuilder::new()))?;
    let graph = graph.unwrap();
    let [a, b, c] = series
        .into_coefficients()
        .map(|root| graph.expression(root));
    Ok(SharedLaurentSeries {
        coefficients: [a?, b?, c?],
    })
}

/// Select branches of a previously generated complete shared expression.
///
/// Probes affect conditions only. All surviving bindings are substituted, so
/// the result is a regular native Symbolica Atom with parametric branch bodies.
/// Undecided conditions remain; a very large result can still exceed `options`.
pub fn select_branch_shared(
    expression: &AliasedAtom,
    rules: &[Replacement],
    options: ExpressionOptions,
) -> Result<Atom, String> {
    if options.max_nodes == 0 || options.max_depth == 0 {
        return Err("expression max_nodes and max_depth must be positive".into());
    }
    let map = FunctionMap::new();
    let mut expansion = Expansion::new(&map, options);
    expansion.aliases = expression
        .get_aliases()
        .iter()
        .map(|(a, b)| (a.clone(), b.clone()))
        .collect();
    expansion.branch_rules = Some(
        rules
            .iter()
            .map(crate::branch_selection::exact_replacement)
            .collect(),
    );
    // Compound patterns can match across binding boundaries. Preserve their
    // native semantics by probing fully expanded predicates instead. The
    // current upstream API exposes no matching-settings getters/PartialEq;
    // compare its derived Debug representation conservatively against default.
    // Python uses an explicit True condition and rhs_cache_size=100 by default;
    // these do not change matching scope and must not disable DAG sampling.
    let ordinary = format!("{:?}", symbolica::id::MatchSettings::default());
    if rules.iter().all(|r| {
        matches!(&r.pat, symbolica::id::Pattern::Literal(a) if matches!(a.as_view(), AtomView::Var(_)))
            && matches!(r.conditions, None | Some(symbolica::id::Condition::True))
            && format!("{:?}", r.match_settings.clone().rhs_cache_size(0)) == ordinary
    }) {
        expansion.graph_probe = Some(NumericProbe::new());
    }
    Ok(expansion
        .expand(expression.get_root().as_view(), 0, 0)?
        .atom
        .as_ref()
        .clone())
}

/// Evaluate a predicate through the stored DAG, without unfolding other arms.
/// Only temporary samples receive replacements. All arithmetic stays in native
/// Symbolica arbitrary precision; a precision increase invalidates every cached
/// sample so no lower-precision intermediate can leak into a higher-precision
/// decision.
pub(super) struct NumericProbe {
    bits: u32,
    values: HashMap<Atom, Complex<Float>>,
    memo: HashMap<(Atom, bool), Option<Complex<Float>>>,
    active: HashSet<(Atom, bool)>,
}

impl NumericProbe {
    fn new() -> Self {
        Self {
            bits: 256,
            values: HashMap::new(),
            memo: HashMap::new(),
            active: HashSet::new(),
        }
    }

    pub(super) fn truth(
        &mut self,
        condition: AtomView<'_>,
        aliases: &HashMap<Atom, Atom>,
        rules: &[Replacement],
    ) -> Option<bool> {
        loop {
            match self.evaluate(condition, aliases, rules, 0, true) {
                Ok(value) => return value.map(|v| !v.is_zero()),
                Err(bits) => {
                    self.bits = bits;
                    self.values.clear();
                    self.memo.clear();
                    self.active.clear();
                }
            }
        }
    }

    fn evaluate(
        &mut self,
        expression: AtomView<'_>,
        aliases: &HashMap<Atom, Atom>,
        rules: &[Replacement],
        depth: usize,
        apply: bool,
    ) -> Result<Option<Complex<Float>>, u32> {
        if depth > 4096 {
            return Ok(None);
        }
        let key = (expression.to_owned(), apply);
        if let Some(value) = self.memo.get(&key) {
            return Ok(value.clone());
        }
        if !self.active.insert(key.clone()) {
            return Ok(None);
        }
        let value = if let Some(body) = aliases.get(&key.0) {
            self.evaluate(body.as_view(), aliases, rules, depth + 1, true)?
        } else {
            let sample = if apply {
                expression.replace_multiple(rules)
            } else {
                expression.to_owned()
            };
            if let Some(bits) = crate::branch_selection::numeric_precision(sample.as_view())
                && bits > self.bits
            {
                return Err(bits);
            }
            if let AtomView::Fun(f) = sample.as_view()
                && f.get_symbol() == Symbol::IF
                && f.get_nargs() == 3
            {
                let args = f.iter().collect::<Vec<_>>();
                if let Some(condition) = self.evaluate(args[0], aliases, rules, depth + 1, false)? {
                    self.evaluate(
                        args[if condition.is_zero() { 2 } else { 1 }],
                        aliases,
                        rules,
                        depth + 1,
                        false,
                    )?
                } else {
                    None
                }
            } else {
                let mut references = Vec::new();
                sample.visitor(&mut |v| {
                    if aliases.contains_key(&v.to_owned()) {
                        references.push(v.to_owned());
                        false
                    } else {
                        true
                    }
                });
                let mut available = true;
                for reference in references {
                    // The probe rules apply to the *definitions* of aliases as
                    // well, but never to the returned symbolic body.
                    if self
                        .evaluate(reference.as_view(), aliases, rules, depth + 1, false)?
                        .is_none()
                    {
                        available = false;
                        break;
                    }
                }
                if available {
                    sample
                        .evaluate_with_prec(&self.values, self.bits)
                        .ok()
                        .filter(|v: &Complex<Float>| v.is_finite())
                } else {
                    None
                }
            }
        };
        self.active.remove(&key);
        if let Some(value) = &value
            && aliases.contains_key(&key.0)
        {
            self.values.insert(key.0.clone(), value.clone());
        }
        self.memo.insert(key, value.clone());
        Ok(value)
    }
}

pub(super) struct GraphBuilder {
    pub(super) id: usize,
    bodies: Vec<(Atom, Atom)>,
    interned: HashMap<Atom, Atom>,
}

impl GraphBuilder {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        Self {
            id: NEXT.fetch_add(1, Ordering::Relaxed),
            bodies: Vec::new(),
            interned: HashMap::new(),
        }
    }

    pub(super) fn intern(&mut self, body: Atom) -> Atom {
        if matches!(body.as_view(), AtomView::Num(_) | AtomView::Var(_)) {
            return body;
        }
        if let Some(alias) = self.interned.get(&body) {
            return alias.clone();
        }
        // Only propagate attributes which Symbolica itself proves for the body.
        let symbol = if body.is_positive().is_true() {
            alias_symbols()[2]
        } else if body.is_real().is_true() {
            alias_symbols()[1]
        } else {
            alias_symbols()[0]
        };
        let alias = symbol.call((self.id, self.bodies.len()));
        self.interned.insert(body.clone(), alias.clone());
        self.bodies.push((alias.clone(), body));
        alias
    }

    fn expression(&self, root: Atom) -> Result<AliasedAtom, String> {
        let mut result = AliasedAtom::from(root);
        // Retain only the bindings reachable from this coefficient. Each body
        // still includes both arms of every unresolved conditional.
        let mut reached = HashSet::new();
        let mut pending = vec![result.get_root().clone()];
        let definitions: HashMap<_, _> =
            self.bodies.iter().map(|(a, b)| (a.as_view(), b)).collect();
        while let Some(atom) = pending.pop() {
            atom.visitor(&mut |view| {
                if let Some(body) = definitions.get(&view) {
                    if reached.insert(view.to_owned()) {
                        pending.push((*body).clone());
                    }
                    false
                } else {
                    true
                }
            });
        }
        for (alias, body) in &self.bodies {
            if reached.contains(alias) {
                let mut unresolved = None;
                body.visitor(&mut |v| {
                    if let Some(s) = v.get_symbol()
                        && (s.get_name().contains("::__olo_")
                            || matches!(
                                s.get_name(),
                                "oneloopmaster::A0"
                                    | "oneloopmaster::B0"
                                    | "oneloopmaster::dB0"
                                    | "oneloopmaster::C0"
                                    | "oneloopmaster::D0"
                            ))
                    {
                        unresolved = Some(s);
                    }
                    unresolved.is_none()
                });
                if let Some(symbol) = unresolved {
                    return Err(format!("incomplete shared definition contains {symbol}"));
                }
                result.register_alias(alias.clone(), body.clone());
            }
        }
        Ok(result)
    }
}
