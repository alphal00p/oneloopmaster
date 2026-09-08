//! Exact Symbolica expressions for one-loop integrals.
//!
//! All masses and renormalization scales passed to this crate are squared
//! quantities. Laurent coefficients are stored in the OneLOop order
//! `[epsilon^0, epsilon^-1, epsilon^-2]`.
#![forbid(unsafe_code)]

mod evaluators;
mod expressions;
mod initialization;
mod sheet_exact;
mod triangle_hv;
pub use evaluators::{
    JitEvaluator, ScalarEvaluator, evaluate, evaluate_batch, jit_settings, rebuild_cached_evaluator,
};
pub use expressions::OneLoopExpressions;
pub use initialization::{initialize, is_initialized};

symbolica::initialize!(
    initialization::from_symbolica,
    "symbolica::special_functions"
);
mod box_complex;
mod box_integral;
mod masters;
mod triangle;
mod two_point;

pub use box_integral::*;
use masters::register_masters;
pub use masters::{A0, B0, C0, D0, ScalarIntegral, dB0};
pub use triangle::*;
pub use two_point::*;

use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    coefficient::Coefficient,
    evaluate::FunctionMap,
    transcendental::TranscendentalFunctions,
};

use sheet_exact::{SheetAtom, divided_difference as sheet_dilog_divided_difference};
use sheet_exact::{SheetAtom as ExactSheetAtom, divided_difference as exact_sheet_difference};

fn native_function_options() -> symbolica::evaluate::FunctionRegistrationOptions {
    symbolica::evaluate::FunctionRegistrationOptions::new()
        .inlining(symbolica::evaluate::InliningPolicy::Never)
}

/// The finite, simple-pole, and double-pole coefficients of an integral.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LaurentSeries {
    coefficients: [Atom; 3],
}

/// A compact Laurent series plus the transparent native definitions needed to evaluate it.
#[derive(Clone, Debug)]
pub struct MappedLaurentSeries {
    series: LaurentSeries,
    function_map: FunctionMap,
}

impl MappedLaurentSeries {
    /// Returns the compact coefficients in `[epsilon^0, epsilon^-1, epsilon^-2]` order.
    pub fn coefficients(&self) -> &[Atom; 3] {
        self.series.coefficients()
    }

    /// Returns the Symbolica function definitions referenced by the coefficients.
    pub fn function_map(&self) -> &FunctionMap {
        &self.function_map
    }

    /// Consumes the result into its compact series and transparent function map.
    pub fn into_parts(self) -> (LaurentSeries, FunctionMap) {
        (self.series, self.function_map)
    }
}

impl LaurentSeries {
    fn new(finite: Atom, simple_pole: Atom, double_pole: Atom) -> Self {
        Self {
            coefficients: [finite, simple_pole, double_pole],
        }
    }

    /// Returns the coefficients in the order `[epsilon^0, epsilon^-1, epsilon^-2]`.
    pub fn coefficients(&self) -> &[Atom; 3] {
        &self.coefficients
    }

    /// Consumes the series and returns its coefficients.
    pub fn into_coefficients(self) -> [Atom; 3] {
        self.coefficients
    }
}

/// Error returned for tadpole tensor ranks above four.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct UnsupportedRank(pub usize);

fn if_nonzero(condition: &Atom, value: Atom) -> Atom {
    Symbol::IF.call((condition, value, 0))
}

fn if_nonzero_else(condition: &Atom, nonzero: Atom, zero: Atom) -> Atom {
    Symbol::IF.call((condition, nonzero, zero))
}

// Test exact zeros without squaring small values or expanding complex norms.
fn any_nonzero(conditions: impl IntoIterator<Item = Atom>) -> Atom {
    conditions
        .into_iter()
        .fold(Atom::num(0), |predicate, condition| {
            if_nonzero_else(&condition, Atom::num(1), predicate)
        })
}

/// Preserve an arithmetic subexpression when forming exact axis predicates.
///
/// Native conjugation distributes over sums and products during normalization.
/// Comparing that expanded result with a separately rounded value can turn an
/// exact real-axis test into a spurious nonzero. The native identity if(z,z,0)
/// keeps both uses grouped, without a callback or a numerical tolerance.
fn grouped(value: &Atom) -> Atom {
    Symbol::IF.call((value, value, 0))
}

fn off_real_axis(value: &Atom) -> Atom {
    let value = grouped(value);
    &value - value.conj()
}

/// Principal logarithm continued to the lower lip of the negative real axis.
///
/// This is the Feynman branch selected by squared masses with non-positive
/// imaginary parts. It is built entirely from Symbolica's native operations.
fn physical_log(value: &Atom) -> Atom {
    // Complex conjugation flips IEEE signed zero on the real axis, so the
    // nested-conjugation identity alone cannot prescribe its branch lip. Select
    // the lower lip explicitly for real arguments and retain the principal
    // logarithm away from the axis.  The IF is native and lazy, so the
    // inactive branch may contain singular expressions without affecting
    // evaluation of the active coefficient.
    let off_axis = off_real_axis(value);
    let value = grouped(value);
    let absolute = Symbol::ABS.call((&value,));
    let negative = Symbol::IF.call((&value - &absolute, Atom::num(1), Atom::num(0)));
    let lower_lip = absolute.log() - Atom::num(-1).sqrt() * Symbol::PI.to_atom() * negative;
    Symbol::IF.call((&off_axis, value.log(), lower_lip))
}

/// Principal logarithm with an explicitly upper real-axis lip. IEEE signed
/// zero is not a reliable substitute for the analytic root's prescribed lip.
fn upper_log(value: &Atom) -> Atom {
    let off_axis = off_real_axis(value);
    let value = grouped(value);
    let absolute = Symbol::ABS.call((&value,));
    let negative = if_nonzero_else(&(&value - &absolute), Atom::num(1), Atom::num(0));
    let upper_lip = absolute.log() + sheet_exact::i() * Symbol::PI.to_atom() * negative;
    if_nonzero_else(&off_axis, value.log(), upper_lip)
}

fn physical_dilog(value: &Atom) -> Atom {
    value.conj().polylog(2).conj()
}

fn physical_sqrt(value: &Atom) -> Atom {
    let off_axis = off_real_axis(value);
    let value = grouped(value);
    let absolute = Symbol::ABS.call((&value,));
    let real_axis = if_nonzero_else(
        &(&value - &absolute),
        -sheet_exact::i() * absolute.sqrt(),
        absolute.sqrt(),
    );
    if_nonzero_else(&off_axis, value.sqrt(), real_axis)
}

fn log_over_one_minus(value: &Atom) -> Atom {
    let difference = 1 - value;
    if_nonzero_else(
        &difference,
        physical_log(value) / &difference,
        Atom::num(-1),
    )
}

fn dilog_complement(value: &Atom) -> Atom {
    physical_dilog(&(1 - value))
}

fn dilog_divided_difference(first: &Atom, second: &Atom) -> Atom {
    let difference = first - second;
    let unequal = (dilog_complement(first) - dilog_complement(second)) / &difference;
    if_nonzero_else(&difference, unequal, log_over_one_minus(first))
}

fn x_log_ratio(value: &Atom, mu_squared: &Atom) -> Atom {
    if value.is_zero() {
        Atom::num(0)
    } else {
        if_nonzero(value, value * physical_log(&(value / mu_squared)))
    }
}

#[cfg(test)]
mod tests;
