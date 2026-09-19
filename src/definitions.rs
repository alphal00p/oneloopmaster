//! Keep the exact definitions we register, independently of Symbolica internals.
use std::{collections::HashMap, sync::OnceLock};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    evaluate::{EvaluationError, FunctionMap as SymbolicaFunctionMap, FunctionRegistrationOptions},
};

#[path = "jit_primitives.rs"]
mod jit_primitives;

type Definition = (Vec<Indeterminate>, Atom, FunctionRegistrationOptions);

#[derive(Clone, Debug, Default)]
pub(crate) struct FunctionMap {
    inner: SymbolicaFunctionMap,
    definitions: HashMap<(Symbol, Vec<Atom>), Definition>,
    compiled: OnceLock<(SymbolicaFunctionMap, bool)>,
    jit_compiled: OnceLock<(SymbolicaFunctionMap, bool)>,
    tags: HashMap<Symbol, usize>,
}

impl FunctionMap {
    pub fn new() -> Self {
        jit_primitives::register();
        Self::default()
    }

    // Reserve aliases before cache decoding can introduce them in traversal order.
    pub fn prepare_symbols(&self) {
        let mut symbols = self.tags.keys().copied().collect::<Vec<_>>();
        symbols.sort();
        for symbol in symbols {
            let _ = symbolica::symbol!(format!("{}__explicit_pi", symbol.get_name()));
        }
    }

    pub fn as_symbolica(&self) -> &SymbolicaFunctionMap {
        &self.compiled.get_or_init(|| self.with_explicit_pi(false)).0
    }

    pub fn as_jit_symbolica(&self) -> &SymbolicaFunctionMap {
        &self
            .jit_compiled
            .get_or_init(|| self.with_explicit_pi(true))
            .0
    }

    // Forward pi and normalize closed helper scopes together, avoiding repeated
    // capture-pruning compilation. Upstream now fixes constant lifting indices;
    // these wrappers retain the shared-scope optimization and original ordering.
    // Public signatures and inspected bodies keep their original definitions.
    fn with_explicit_pi(&self, jit: bool) -> (SymbolicaFunctionMap, bool) {
        let primitives = jit.then(jit_primitives::register);
        let pi = Symbol::PI.to_atom();
        // Atom canonical ordering uses symbol IDs. Preserve the helpers' relative
        // order when creating aliases, independently of HashMap iteration order.
        let mut symbols = self.tags.keys().copied().collect::<Vec<_>>();
        symbols.sort();
        let implementations: HashMap<Symbol, Symbol> = symbols
            .into_iter()
            .map(|symbol| {
                (
                    symbol,
                    symbolica::symbol!(format!("{}__explicit_pi", symbol.get_name())),
                )
            })
            .collect();
        // Closed helpers can share formal slots. Padding their signatures keeps
        // nested scopes identical, so upstream need not recursively rebuild each
        // callee to remove unused lexical captures. Keep ordinary lexical scopes
        // for maps that contain free variables.
        let closed = self.definitions.values().all(|(args, body, _)| {
            let mut closed = true;
            body.visitor(&mut |view| {
                if args.iter().any(|arg| arg.as_view() == view) {
                    return false;
                }
                if matches!(view, AtomView::Var(_))
                    && view != pi.as_view()
                    && self.get_definition(view).is_none()
                {
                    closed = false;
                }
                closed
            });
            closed
        });
        let pi_arg = if closed {
            Symbol::PI
        } else {
            symbolica::symbol!("oneloopmaster::__constant_pi")
        };
        let slots = closed.then(|| {
            (0..self
                .definitions
                .values()
                .map(|(a, _, _)| a.len())
                .max()
                .unwrap_or(0))
                .map(|index| symbolica::symbol!(format!("oneloopmaster::__argument_{index}")))
                .collect::<Vec<_>>()
        });
        let mut result = SymbolicaFunctionMap::new();
        for ((symbol, tags), (args, body, options)) in &self.definitions {
            let implementation = implementations[symbol];
            let body = body.replace_map_bottom_up(|view, _, out| {
                if let Some(index) = slots
                    .as_ref()
                    .and_then(|_| args.iter().position(|a| a.as_view() == view))
                {
                    **out = slots.as_ref().unwrap()[index].to_atom();
                } else if view == pi.as_view() {
                    **out = pi_arg.to_atom();
                } else if let Some(callee) = view.get_symbol()
                    && self.get_definition(view).is_some()
                {
                    let mut values = match view {
                        AtomView::Fun(call) => {
                            call.iter().map(|a| a.to_owned()).collect::<Vec<_>>()
                        }
                        _ => Vec::new(),
                    };
                    if let Some(slots) = &slots {
                        values.resize(self.tags[&callee] + slots.len(), Atom::num(0));
                    }
                    values.push(pi_arg.to_atom());
                    **out = implementations[&callee].call(&values);
                } else if let Some(atom) =
                    primitives.and_then(|symbols| jit_primitives::lower(view, symbols))
                {
                    **out = atom;
                }
            });
            let mut explicit_args = slots.as_ref().map_or_else(
                || args.clone(),
                |slots| slots.iter().copied().map(Into::into).collect(),
            );
            explicit_args.push(pi_arg.into());
            result
                .add_tagged_function_with_options(
                    implementation,
                    tags.clone(),
                    explicit_args,
                    body,
                    options.clone(),
                )
                .unwrap();
            let mut values = tags.clone();
            values.extend(args.iter().map(|arg| arg.as_view().to_owned()));
            if let Some(slots) = &slots {
                values.resize(tags.len() + slots.len(), Atom::num(0));
            }
            values.push(Symbol::PI.to_atom());
            // Default inlining resolves this wrapper in the caller's scope.
            result
                .add_tagged_function(
                    *symbol,
                    tags.clone(),
                    args.clone(),
                    implementation.call(&values),
                )
                .unwrap();
        }
        (result, closed)
    }

    pub fn normalize_inputs(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
        jit: bool,
    ) -> (Vec<Atom>, Vec<Atom>) {
        let prepared = if jit {
            &self.jit_compiled
        } else {
            &self.compiled
        };
        if !prepared.get_or_init(|| self.with_explicit_pi(jit)).1 {
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
            return (expressions, parameters.to_vec());
        }
        let primitives = jit.then(jit_primitives::register);
        let normalized: Vec<Atom> = parameters
            .iter()
            .enumerate()
            .map(|(index, parameter)| {
                if parameter == &Symbol::PI.to_atom() {
                    parameter.clone()
                } else {
                    symbolica::symbol!(format!("oneloopmaster::__argument_{index}")).to_atom()
                }
            })
            .collect();
        let arity = self
            .definitions
            .values()
            .map(|(args, _, _)| args.len())
            .max()
            .unwrap_or(0);
        let expressions = expressions
            .iter()
            .map(|expression| {
                expression.replace_map_bottom_up(|view, _, out| {
                    if let Some(index) = parameters.iter().position(|p| p.as_view() == view) {
                        **out = normalized[index].clone();
                    } else if let Some((tags, args, _)) = self.get_definition(view) {
                        let mut values = match view {
                            AtomView::Fun(call) => {
                                call.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>()
                            }
                            _ => Vec::new(),
                        };
                        if values.len() == tags + args.len() {
                            values.resize(tags + arity, Atom::num(0));
                            values.push(Symbol::PI.to_atom());
                            **out = symbolica::symbol!(format!(
                                "{}__explicit_pi",
                                view.get_symbol().unwrap().get_name()
                            ))
                            .call(&values);
                        }
                    } else if let Some(atom) =
                        primitives.and_then(|symbols| jit_primitives::lower(view, symbols))
                    {
                        **out = atom;
                    }
                })
            })
            .collect();
        (expressions, normalized)
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
        self.compiled.take();
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
}

impl From<FunctionMap> for SymbolicaFunctionMap {
    fn from(map: FunctionMap) -> Self {
        map.as_symbolica().clone()
    }
}
