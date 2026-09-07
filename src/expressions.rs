//! Reusable definitions and evaluator construction for mixed integral expressions.
use super::*;
use std::sync::OnceLock;
use symbolica::{
    domains::{float::Complex, rational::Rational},
    evaluate::{EvaluationError, ExpressionEvaluator},
};

/// Native definitions shared by any number of triangle and box expressions.
/// Construct this once when building an amplitude, then pass its function map
/// to Symbolica or use [Self::evaluator] to compile a combination of integrals.
#[derive(Clone, Debug)]
pub struct OneLoopExpressions {
    map: FunctionMap,
}

impl Default for OneLoopExpressions {
    fn default() -> Self {
        Self::new()
    }
}

impl OneLoopExpressions {
    /// Builds the native definitions once per process and reuses them thereafter.
    pub fn new() -> Self {
        static DEFINITIONS: OnceLock<FunctionMap> = OnceLock::new();
        Self {
            map: DEFINITIONS
                .get_or_init(|| {
                    let mut map = FunctionMap::new();
                    sheet_exact::register(&mut map);
                    triangle_hv::register(&mut map);
                    register_triangles(&mut map);
                    register_boxes(&mut map);
                    map
                })
                .clone(),
        }
    }

    /// Returns the shared transparent Symbolica definitions.
    pub fn function_map(&self) -> &FunctionMap {
        &self.map
    }
    /// Transfers the definitions to a Symbolica evaluator builder.
    pub fn into_function_map(self) -> FunctionMap {
        self.map
    }

    /// Constructs a triangle using this context's definitions.
    /// The same unsupported degenerate limits as [crate::c0] apply.
    pub fn c0(&self, p: [&Atom; 3], m: [&Atom; 3], mu: &Atom) -> LaurentSeries {
        let args = p
            .into_iter()
            .chain(m)
            .chain([mu])
            .cloned()
            .collect::<Vec<_>>();
        let branches = core::array::from_fn(|sector| series_call("c0", sector, &args));
        with_triangle_normalization(select_three_series(m, branches))
    }

    /// Constructs a box using this context's definitions.
    /// The same coverage limitations as [crate::d0] apply.
    pub fn d0(&self, p: [&Atom; 6], m: [&Atom; 4], mu: &Atom) -> LaurentSeries {
        let args = p
            .into_iter()
            .chain(m)
            .chain([mu])
            .cloned()
            .collect::<Vec<_>>();
        select_four_series(
            m,
            core::array::from_fn(|sector| series_call("d0", sector, &args)),
        )
    }

    /// Compiles combinations of integral coefficients. The returned evaluator
    /// retains exact rational constants. See the README for stack requirements.
    pub fn evaluator(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
        compile(expressions, parameters, self.map.clone())
    }
}

fn series_call(family: &str, sector: usize, args: &[Atom]) -> LaurentSeries {
    let call = |coefficient| {
        symbolica::symbol!(format!("__olo_{family}_sector_{sector}_{coefficient}")).call(args)
    };
    LaurentSeries::new(call(0), call(1), call(2))
}

fn compile(
    expressions: &[Atom],
    parameters: &[Atom],
    map: FunctionMap,
) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
    let views = expressions.iter().map(Atom::as_view).collect::<Vec<_>>();
    Atom::evaluator_multiple(&views, parameters)
        .function_map(map)
        .direct_translation(true)
        .build()
}

impl MappedLaurentSeries {
    /// Compiles the three coefficients using the included definitions.
    pub fn evaluator(
        &self,
        parameters: &[Atom],
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
        compile(self.coefficients(), parameters, self.function_map.clone())
    }
}
