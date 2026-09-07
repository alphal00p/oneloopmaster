//! Exact Symbolica expressions for one-loop integrals.
//!
//! All masses and renormalization scales passed to this crate are squared
//! quantities. Laurent coefficients are stored in the OneLOop order
//! `[epsilon^0, epsilon^-1, epsilon^-2]`.
#![forbid(unsafe_code)]

mod expressions;
mod sheet_exact;
mod triangle_hv;
pub use expressions::OneLoopExpressions;
mod box_integral;
mod triangle;
mod two_point;

pub use box_integral::*;
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

/// Principal logarithm continued to the lower lip of the negative real axis.
///
/// This is the Feynman branch selected by squared masses with non-positive
/// imaginary parts. It is built entirely from Symbolica's native operations.
fn physical_log(value: &Atom) -> Atom {
    value.conj().log().conj()
}

fn physical_dilog(value: &Atom) -> Atom {
    value.conj().polylog(2).conj()
}

fn physical_sqrt(value: &Atom) -> Atom {
    value.conj().sqrt().conj()
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
