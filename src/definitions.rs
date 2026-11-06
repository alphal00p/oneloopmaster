//! Keep the exact definitions we register, independently of Symbolica internals.
use std::{
    collections::{HashMap, HashSet},
    sync::OnceLock,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    evaluate::{
        EvaluationError, FunctionMap as SymbolicaFunctionMap, FunctionRegistrationOptions,
        InliningPolicy,
    },
};

#[path = "jit_primitives.rs"]
mod jit_primitives;

pub(crate) fn register_callbacks() {
    jit_primitives::register();
}

/// Opaque during Atom normalization, transparent during instruction lowering.
/// Direct translation evaluates the argument into the caller's scope before
/// inlining this identity; conjugation therefore uses the same rounded value.
fn analytic_grouping_symbol() -> Symbol {
    symbolica::symbol!("__olo_analytic_group")
}

/// Preserve grouping and every Taylor component at a zero primal, without
/// introducing conditional scopes or duplicating the argument's instructions.
pub(crate) fn analytic_grouping(expression: &Atom) -> Atom {
    expression.replace_map_bottom_up(|view, _, out| {
        let AtomView::Fun(fun) = view else {
            return;
        };
        if fun.get_symbol() != Symbol::IF || fun.get_nargs() != 3 {
            return;
        }
        let mut args = fun.iter();
        let condition = args.next().unwrap();
        let nonzero = args.next().unwrap();
        let zero = args.next().unwrap();
        if condition == nonzero && zero.is_zero() {
            **out = analytic_grouping_symbol().call((condition,));
        }
    })
}

type Definition = (Vec<Indeterminate>, Atom, FunctionRegistrationOptions);

#[derive(Clone, Debug, Default)]
pub(crate) struct FunctionMap {
    inner: SymbolicaFunctionMap,
    definitions: HashMap<(Symbol, Vec<Atom>), Definition>,
    jit_compiled: OnceLock<SymbolicaFunctionMap>,
    tags: HashMap<Symbol, usize>,
}

impl FunctionMap {
    pub fn new() -> Self {
        jit_primitives::register();
        Self::default()
    }

    pub fn as_symbolica(&self) -> &SymbolicaFunctionMap {
        &self.inner
    }

    pub fn as_jit_symbolica(&self) -> &SymbolicaFunctionMap {
        self.jit_compiled.get_or_init(|| {
            let primitives = jit_primitives::register();
            let mut result = SymbolicaFunctionMap::new();
            for ((symbol, tags), (args, body, options)) in &self.definitions {
                let body = body.replace_map_bottom_up(|view, _, out| {
                    if let Some(atom) = jit_primitives::lower(view, primitives) {
                        **out = atom;
                    }
                });
                result
                    .add_tagged_function_with_options(
                        *symbol,
                        tags.clone(),
                        args.clone(),
                        body,
                        options.clone(),
                    )
                    .unwrap();
            }
            result
        })
    }

    // Upstream registers non-inlined functions and their constants globally.
    // Keep original signatures and symbol ordering; shared argument slots and
    // explicit-pi forwarding are no longer needed to avoid recursive rebuilding.
    pub fn normalize_inputs(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
        jit: bool,
    ) -> (Vec<Atom>, Vec<Atom>) {
        let expressions = if jit {
            let primitives = jit_primitives::register();
            expressions
                .iter()
                .map(|expression| {
                    expression.replace_map_bottom_up(|view, _, out| {
                        if let Some(atom) = jit_primitives::lower(view, primitives) {
                            **out = atom;
                        }
                    })
                })
                .collect()
        } else {
            expressions.to_vec()
        };
        (expressions, parameters.to_vec())
    }

    pub fn add_function<S: Into<Indeterminate>, A: Into<Indeterminate>>(
        &mut self,
        name: S,
        args: Vec<A>,
        body: Atom,
    ) -> Result<(), EvaluationError> {
        self.add_function_with_options(name, args, body, FunctionRegistrationOptions::default())
    }

    pub fn add_function_with_options<S: Into<Indeterminate>, A: Into<Indeterminate>>(
        &mut self,
        name: S,
        args: Vec<A>,
        body: Atom,
        options: FunctionRegistrationOptions,
    ) -> Result<(), EvaluationError> {
        let (symbol, tags) = match name.into() {
            Indeterminate::Symbol(symbol, _) => (symbol, Vec::new()),
            Indeterminate::Function(symbol, call) => (
                symbol,
                call.as_fun_view()
                    .unwrap()
                    .iter()
                    .map(|a| a.to_owned())
                    .collect(),
            ),
        };
        self.add_tagged_function_with_options(symbol, tags, args, body, options)
    }

    pub fn add_tagged_function_with_options<A: Into<Indeterminate>>(
        &mut self,
        symbol: Symbol,
        tags: Vec<Atom>,
        args: Vec<A>,
        body: Atom,
        options: FunctionRegistrationOptions,
    ) -> Result<(), EvaluationError> {
        let args: Vec<Indeterminate> = args.into_iter().map(Into::into).collect();
        // Let Symbolica validate the registration before updating our inspection data.
        self.inner.add_tagged_function_with_options(
            symbol,
            tags.clone(),
            args.clone(),
            body.clone(),
            options.clone(),
        )?;
        self.tags.insert(symbol, tags.len());
        self.definitions
            .insert((symbol, tags), (args, body, options));
        self.jit_compiled.take();
        Ok(())
    }

    pub fn get_definition(&self, call: AtomView<'_>) -> Option<(usize, &[Indeterminate], &Atom)> {
        let symbol = call.get_symbol()?;
        let count = *self.tags.get(&symbol)?;
        let tags = match call {
            AtomView::Var(_) if count == 0 => Vec::new(),
            AtomView::Fun(f) if f.get_nargs() >= count => {
                f.iter().take(count).map(|a| a.to_owned()).collect()
            }
            _ => return None,
        };
        self.definitions
            .get(&(symbol, tags))
            .map(|(args, body, _)| (count, args.as_slice(), body))
    }

    /// Byte size of every registered body, largest first.
    pub(crate) fn body_sizes(&self) -> Vec<(String, usize)> {
        let mut sizes: Vec<(String, usize)> = self
            .definitions
            .iter()
            .map(|((symbol, _), (_, body, _))| {
                (
                    symbol.get_name().to_string(),
                    body.as_view().get_byte_size(),
                )
            })
            .collect();
        sizes.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        sizes
    }

    /// Prepare algebraic grouping identities for arbitrary Dual directions.
    pub(crate) fn with_analytic_grouping(&self) -> Result<Self, String> {
        let mut map = Self::new();
        let registered: HashSet<Symbol> =
            self.definitions.keys().map(|(symbol, _)| *symbol).collect();
        for ((symbol, tags), (args, body, options)) in &self.definitions {
            // Retained scalar helpers become one helper per Taylor component.
            // Inline bounded leaves and pure forwarding aliases before Dualizer,
            // so Symbolica can share/prune their arithmetic in the caller. Keep
            // compound call graphs outlined to bound preparation and IR growth.
            let leaf = body.as_view().get_byte_size() <= 256
                && body
                    .get_all_symbols(true)
                    .iter()
                    .all(|symbol| !registered.contains(symbol));
            let forwarding = match body.as_view() {
                AtomView::Fun(fun) if registered.contains(&fun.get_symbol()) => fun
                    .iter()
                    .all(|argument| matches!(argument, AtomView::Var(_) | AtomView::Num(_))),
                _ => false,
            };
            let options = if leaf || forwarding {
                options.clone().inlining(InliningPolicy::Always)
            } else {
                options.clone()
            };
            map.add_tagged_function_with_options(
                *symbol,
                tags.clone(),
                args.clone(),
                analytic_grouping(body),
                options,
            )
            .map_err(|error| error.to_string())?;
        }
        let group = analytic_grouping_symbol();
        if !map.definitions.contains_key(&(group, Vec::new())) {
            let argument = symbolica::symbol!("__olo_analytic_group_argument");
            map.add_function_with_options(
                group,
                vec![argument],
                argument.to_atom(),
                FunctionRegistrationOptions::new()
                    .inlining(symbolica::evaluate::InliningPolicy::Always),
            )
            .map_err(|error| error.to_string())?;
        }
        Ok(map)
    }

    /// Rebuild the recorded map with one untagged helper's body replaced.
    pub(crate) fn with_replaced_body(&self, name: Symbol, body: Atom) -> Result<Self, String> {
        if !self.definitions.contains_key(&(name, Vec::new())) {
            return Err(format!("missing neighborhood helper {}", name.get_name()));
        }
        let mut map = Self::new();
        for ((symbol, tags), (args, original, options)) in &self.definitions {
            let replacement = if *symbol == name && tags.is_empty() {
                &body
            } else {
                original
            };
            map.add_tagged_function_with_options(
                *symbol,
                tags.clone(),
                args.clone(),
                replacement.clone(),
                options.clone(),
            )
            .map_err(|error| error.to_string())?;
        }
        Ok(map)
    }
}

impl From<FunctionMap> for SymbolicaFunctionMap {
    fn from(map: FunctionMap) -> Self {
        map.as_symbolica().clone()
    }
}
