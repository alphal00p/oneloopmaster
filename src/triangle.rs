//! Triangle sectors and their native expressions.
use super::*;

pub(crate) fn negative_quadratic_roots(a: &Atom, b: &Atom, c: &Atom) -> (Atom, Atom) {
    let discriminant = physical_sqrt(&(b * b - Atom::num(4) * a * c));
    (
        (b + &discriminant) / (Atom::num(2) * a),
        (b - discriminant) / (Atom::num(2) * a),
    )
}

pub(crate) fn r_function(value: &Atom) -> (Atom, Atom) {
    // OneLOop selects the root outside the unit circle for non-real q and
    // real |q|>2, but the lower-half-plane root on the real interval [-2,2].
    // Choosing the reciprocal preserves q=r+1/r but changes continued Li2 terms.
    let difference = (value - 2).sqrt() * (value + 2).sqrt();
    let outside = (value + &difference) / 2;
    let discriminant: Atom = value * value - 4;
    // Keep the real-axis formulas on their exact axes: multiplying two
    // numerically evaluated negative-real square roots can leak a tiny
    // imaginary component and incorrectly select a continuation sheet.
    let real_inside = (value - sheet_exact::i() * (-&discriminant).sqrt()) / 2;
    let real_outside = (value + sheet_exact::sign_nonnegative(value) * discriminant.sqrt()) / 2;
    let root = if_nonzero_else(
        &sheet_exact::im(value),
        outside,
        if_nonzero_else(
            &sheet_exact::negative(&discriminant),
            real_inside,
            real_outside,
        ),
    );
    let inverse = Atom::num(1) / &root;
    (root.clone(), root - inverse)
}

pub(crate) fn triangle_ir_one(mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let factor = Atom::num(1) / (Atom::num(2) * mass);
    LaurentSeries::new(
        &factor * (2 + physical_log(&(mass / mu_squared))),
        -factor,
        Atom::num(0),
    )
}

pub(crate) fn triangle_ir_two(momentum: &Atom, mass: &Atom, mu_squared: &Atom) -> LaurentSeries {
    let difference = mass - momentum;
    let q_mass = SheetAtom::lower(mass.clone());
    let q_difference = SheetAtom::lower(difference.clone());
    let log_mass = q_mass.scale_positive(&(Atom::num(1) / mu_squared)).log();
    let ratio = &q_mass / &q_difference;
    let logarithm = ratio.log();
    let double_pole = Atom::num(1) / 2;
    let single_pole = &logarithm - &double_pole * &log_mass;
    let finite_seed = Atom::num(1) / 24 * Symbol::PI.to_atom() * Symbol::PI.to_atom()
        + &logarithm * &logarithm / 2
        - ratio.dilog();
    let factor = -Atom::num(1) / difference;
    LaurentSeries::new(
        &factor * (finite_seed + (&double_pole * &log_mass / 2 - logarithm) * log_mass),
        &factor * single_pole,
        factor * double_pole,
    )
}

pub(crate) fn triangle_ir_three(
    momentum_2: &Atom,
    momentum_3: &Atom,
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let difference_13 = mass - momentum_3;
    let difference_23 = mass - momentum_2;
    let q13 = SheetAtom::lower(difference_13.clone());
    let q23 = SheetAtom::lower(difference_23);
    let q_mass = SheetAtom::lower(mass.clone());
    let x1 = &q23 / &q_mass;
    let x2 = &q13 / &q_mass;
    let single_pole = -(&q23 / &q13).log_over_one_minus() / &difference_13;
    LaurentSeries::new(
        -sheet_dilog_divided_difference(&x1, &x2) / mass
            - &single_pole
                * (x1.log() + x2.log() + q_mass.scale_positive(&(Atom::num(1) / mu_squared)).log()),
        single_pole,
        Atom::num(0),
    )
}

pub(crate) fn triangle_ir_four(
    momentum: &Atom,
    mass_2: &Atom,
    mass_3: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let sqrt_2 = mass_2.sqrt();
    let sqrt_3 = mass_3.sqrt();
    let (r23, _) = r_function(&((mass_2 + mass_3 - momentum) / (&sqrt_2 * &sqrt_3)));
    let q23 = SheetAtom::lower(r23.clone());
    let coefficient = q23.log_over_one_minus() * &r23 / ((Atom::num(1) + &r23) * &sqrt_2 * &sqrt_3);
    let ratio = &SheetAtom::upper(sqrt_3.clone()) / &SheetAtom::upper(sqrt_2.clone());
    let finite = &coefficient * (physical_log(&(mass_3 / mu_squared)) - q23.log())
        - sheet_dilog_divided_difference(&(&ratio * &q23), &(&ratio / &q23)) / mass_2
        + sheet_dilog_divided_difference(&(&q23 * &q23), &SheetAtom::lower(Atom::num(1))) * &r23
            / (&sqrt_2 * &sqrt_3);
    LaurentSeries::new(finite, -coefficient, Atom::num(0))
}

pub(crate) fn triangle_finite_massless(p1: &Atom, p2: &Atom, p3: &Atom) -> Atom {
    let r23 = -p1;
    let r24 = -p3;
    let r34 = -p2;
    let a = &r34 * &r24;
    let b = &r24 + &r34 - &r23;
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &Atom::num(1));
    let qx1 = SheetAtom::with_real_axis_sign(root_1.clone(), true);
    let qx2 = SheetAtom::with_real_axis_sign(root_2.clone(), false);
    let q23 = SheetAtom::lower(r23);
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let x1_34 = &qx1 * &q34;
    let x2_34 = &qx2 * &q34;
    let x1_24 = &qx1 * &q24;
    let x2_24 = &qx2 * &q24;
    let ratio = &qx1 / &qx2;
    let root_product = &qx1 * &qx2;
    let log_ratio = ratio.log_over_one_minus();
    (sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34
        + sheet_dilog_divided_difference(&x1_24, &x2_24) * &r24
        - &log_ratio * root_product.log() / (Atom::num(2) * &root_2)
        - log_ratio * q23.log() / root_2)
        / a
}

pub(crate) fn square_root_magnitude(mass: &Atom) -> Atom {
    (mass * mass.conj()).sqrt().sqrt()
}

pub(crate) fn triangle_finite_one_mass(p1: &Atom, p2: &Atom, p3: &Atom, mass: &Atom) -> Atom {
    let p2_internal = p1;
    let p3_internal = p2;
    let p23 = p3;
    let sqrt_mass = mass.sqrt();
    let sqrt_real = square_root_magnitude(mass);
    let r23 = -p2_internal / (&sqrt_real * &sqrt_real);
    let r24 = (mass - p23) / (&sqrt_real * &sqrt_mass);
    let r34 = (mass - p3_internal) / (&sqrt_real * &sqrt_mass);
    let a = &r34 * &r24 - &r23;
    let b = &r24 / &sqrt_real + &r34 / &sqrt_real - &r23 / &sqrt_mass;
    let c = Atom::num(1) / (&sqrt_real * &sqrt_real);
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &c);
    let qx1 = SheetAtom::with_real_axis_sign(root_1.clone(), true);
    let qx2 = SheetAtom::with_real_axis_sign(root_2.clone(), true);
    let q23 = SheetAtom::lower(r23.clone());
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let qmass = SheetAtom::lower(sqrt_mass.clone());
    let root_product = &qx1 * &qx2;
    let mass_product = &qmass * &qmass;
    let root_over_mass = &root_product / &mass_product;
    let ratio = &qx1 / &qx2;
    let mut finite =
        root_over_mass.log() / 2 + q23.scale_positive(&(&sqrt_real * &sqrt_real)).log();
    let x1_mass = &qx1 * &qmass;
    let x2_mass = &qx2 * &qmass;
    finite = -finite * ratio.log_over_one_minus() / &root_2
        - sheet_dilog_divided_difference(&x1_mass, &x2_mass) * &sqrt_mass;
    let scaled_34 = q34.scale_positive(&sqrt_real);
    let x1_34 = &qx1 * &scaled_34;
    let x2_34 = &qx2 * &scaled_34;
    finite += if_nonzero(
        &r34,
        sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34 * &sqrt_real,
    );
    let scaled_24 = q24.scale_positive(&sqrt_real);
    let x1_24 = &qx1 * &scaled_24;
    let x2_24 = &qx2 * &scaled_24;
    finite += if_nonzero(
        &r24,
        sheet_dilog_divided_difference(&x1_24, &x2_24) * &r24 * &sqrt_real,
    );
    if_nonzero(&a, finite / (&a * &sqrt_real * &sqrt_real * sqrt_mass))
}

pub(crate) fn triangle_finite_two_masses(
    p1: &Atom,
    p2: &Atom,
    p3: &Atom,
    mass_2_input: &Atom,
    mass_3_input: &Atom,
) -> Atom {
    let p2_internal = p3;
    let p3_internal = p1;
    let p23 = p2;
    let mass_2 = mass_3_input;
    let mass_4 = mass_2_input;
    let sqrt_2 = mass_2.sqrt();
    let sqrt_3 = square_root_magnitude(mass_2);
    let sqrt_4 = mass_4.sqrt();
    let r23 = (mass_2 - p2_internal) / (&sqrt_2 * &sqrt_3);
    let k24 = (mass_2 + mass_4 - p23) / (&sqrt_2 * &sqrt_4);
    let r34 = (mass_4 - p3_internal) / (&sqrt_3 * &sqrt_4);
    let (r24, d24) = r_function(&k24);
    let a = &r34 / &r24 - &r23;
    let b = -&d24 / &sqrt_3 + &r34 / &sqrt_2 - &r23 / &sqrt_4;
    let c = (&sqrt_4 / &sqrt_2 - &r24) / (&sqrt_3 * &sqrt_4);
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &c);
    let qx1 = SheetAtom::with_real_axis_sign(root_1.clone(), true);
    let qx2 = SheetAtom::with_real_axis_sign(root_2.clone(), true);
    let q23 = SheetAtom::lower(r23.clone());
    let q24 = SheetAtom::lower(r24.clone());
    let q34 = SheetAtom::lower(r34.clone());
    let qm2 = SheetAtom::lower(sqrt_2.clone());
    let qm3 = SheetAtom::lower(sqrt_3.clone());
    let qm4 = SheetAtom::lower(sqrt_4.clone());
    let y1 = &qx1 / &q24;
    let y2 = &qx2 / &q24;
    let y1m2 = &y1 * &qm2;
    let y2m2 = &y2 * &qm2;
    let mut finite = sheet_dilog_divided_difference(&y1m2, &y2m2) * &sqrt_2 / &r24;
    let y_ratio = &y1 / &y2;
    let y_product = &y1 * &y2;
    let m2_product = &qm2 * &qm2;
    let y_scaled = &y_product / &m2_product;
    let x_ratio = &qx1 / &qx2;
    let x_product = &qx1 * &qx2;
    let m4_product = &qm4 * &qm4;
    let x_scaled = &x_product / &m4_product;
    let root_term = (y_ratio.log_over_one_minus() * y_scaled.log()
        - x_ratio.log_over_one_minus() * x_scaled.log())
        / (Atom::num(2) * &root_2);
    finite += if_nonzero(&root_2, root_term);
    let x1m4 = &qx1 * &qm4;
    let x2m4 = &qx2 * &qm4;
    finite -= sheet_dilog_divided_difference(&x1m4, &x2m4) * &sqrt_4;
    let q23m3 = &q23 * &qm3;
    let scale_23 = &q23m3 / &q24;
    let x1_23 = &qx1 * &scale_23;
    let x2_23 = &qx2 * &scale_23;
    finite -= if_nonzero(
        &r23,
        sheet_dilog_divided_difference(&x1_23, &x2_23) * &r23 * &sqrt_3 / &r24,
    );
    let scale_34 = &q34 * &qm3;
    let x1_34 = &qx1 * &scale_34;
    let x2_34 = &qx2 * &scale_34;
    finite += if_nonzero(
        &r34,
        sheet_dilog_divided_difference(&x1_34, &x2_34) * &r34 * &sqrt_3,
    );
    if_nonzero(&a, finite / (&a * sqrt_2 * sqrt_3 * sqrt_4))
}

pub(crate) fn triangle_finite_three_masses(momenta: [&Atom; 3], masses: [&Atom; 3]) -> Atom {
    let sqrt_1 = masses[0].sqrt();
    let sqrt_2 = masses[1].sqrt();
    let sqrt_3 = masses[2].sqrt();
    let k12 = (masses[0] + masses[1] - momenta[0]) / (&sqrt_1 * &sqrt_2);
    let k13 = (masses[0] + masses[2] - momenta[2]) / (&sqrt_1 * &sqrt_3);
    let k23 = (masses[1] + masses[2] - momenta[1]) / (&sqrt_2 * &sqrt_3);
    let (r12, _) = r_function(&k12);
    let (r13, d13) = r_function(&k13);
    let (r23, _) = r_function(&k23);
    let a = &sqrt_2 / &sqrt_3 - &k23 + &r13 * (&k12 - &sqrt_2 / &sqrt_1);
    let b = &d13 / &sqrt_2 + &k12 / &sqrt_3 - &k23 / &sqrt_1;
    let c = (&sqrt_1 / &sqrt_3 - Atom::num(1) / &r13) / (&sqrt_1 * &sqrt_2);
    let (root_1, root_2) = negative_quadratic_roots(&a, &b, &c);

    use sheet_exact::{SheetAtom as Q, divided_difference as dd};
    let qx1 = Q::upper(root_1.clone());
    let qx2 = Q::upper(root_2.clone());
    let [q12, q13, q23, qm1, qm2, qm3] =
        [&r12, &r13, &r23, &sqrt_1, &sqrt_2, &sqrt_3].map(|r| Q::lower(r.clone()));
    let z1 = &qx1 * &qm2;
    let z2 = &qx2 * &qm2;
    let mut finite = (dd(&(&z1 * &q12), &(&z2 * &q12)) * &r12
        + dd(&(&z1 / &q12), &(&z2 / &q12)) / &r12)
        * &sqrt_2;
    let temp = &q13 * &qm2;
    let t1 = &qx1 * &temp;
    let t2 = &qx2 * &temp;
    finite -= (dd(&(&t1 * &q23), &(&t2 * &q23)) * &r23 + dd(&(&t1 / &q23), &(&t2 / &q23)) / &r23)
        * &r13
        * &sqrt_2;
    let y1 = &qx1 * &q13;
    let y2 = &qx2 * &q13;
    finite += dd(&(&y1 * &qm3), &(&y2 * &qm3)) * &r13 * &sqrt_3
        - dd(&(&qx1 * &qm1), &(&qx2 * &qm1)) * &sqrt_1;
    let root_term = ((&y1 / &y2).log_over_one_minus()
        * (y1.log() + y2.log() - Atom::num(2) * qm3.log())
        - (&qx1 / &qx2).log_over_one_minus() * (qx1.log() + qx2.log() - Atom::num(2) * qm1.log()))
        / (Atom::num(2) * &root_2);
    finite += if_nonzero(&root_2, root_term);
    if_nonzero(&a, finite / (&a * sqrt_1 * sqrt_2 * sqrt_3))
}

pub(crate) fn triangle_massless(momenta: [&Atom; 3], mu_squared: &Atom) -> LaurentSeries {
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    let one_scale = |momentum: &Atom| {
        let logarithm = physical_log(&(-momentum / mu_squared));
        LaurentSeries::new(
            (&logarithm * &logarithm / 2 - &pi_squared / 12) / momentum,
            -&logarithm / momentum,
            Atom::num(1) / momentum,
        )
    };
    let two_scale = |first: &Atom, second: &Atom| {
        let first_log = physical_log(&(-first / mu_squared));
        let second_log = physical_log(&(-second / mu_squared));
        let difference = first - second;
        // A logarithm of the ratio loses the individual Feynman lips when
        // the real invariants have opposite signs.
        let single_pole = if_nonzero_else(
            &difference,
            (&second_log - &first_log) / &difference,
            -Atom::num(1) / first,
        );
        LaurentSeries::new(
            -&single_pole * (second_log + first_log) / 2,
            single_pole,
            Atom::num(0),
        )
    };
    let finite = LaurentSeries::new(
        triangle_finite_massless(momenta[0], momenta[1], momenta[2]),
        Atom::num(0),
        Atom::num(0),
    );
    let zero = LaurentSeries::new(Atom::num(0), Atom::num(0), Atom::num(0));

    let branches = [
        zero,
        one_scale(momenta[2]),
        one_scale(momenta[1]),
        two_scale(momenta[1], momenta[2]),
        one_scale(momenta[0]),
        two_scale(momenta[0], momenta[2]),
        two_scale(momenta[0], momenta[1]),
        finite,
    ];
    select_three_series(momenta, branches)
}

pub(crate) fn select_three_series(
    conditions: [&Atom; 3],
    branches: [LaurentSeries; 8],
) -> LaurentSeries {
    let coefficients = branches.map(LaurentSeries::into_coefficients);
    let select = |index: usize| {
        let with_first = if_nonzero_else(
            conditions[1],
            if_nonzero_else(
                conditions[2],
                coefficients[7][index].clone(),
                coefficients[6][index].clone(),
            ),
            if_nonzero_else(
                conditions[2],
                coefficients[5][index].clone(),
                coefficients[4][index].clone(),
            ),
        );
        let without_first = if_nonzero_else(
            conditions[1],
            if_nonzero_else(
                conditions[2],
                coefficients[3][index].clone(),
                coefficients[2][index].clone(),
            ),
            if_nonzero_else(
                conditions[2],
                coefficients[1][index].clone(),
                coefficients[0][index].clone(),
            ),
        );
        if_nonzero_else(conditions[0], with_first, without_first)
    };
    LaurentSeries::new(select(0), select(1), select(2))
}

pub(crate) fn if_series(
    condition: &Atom,
    nonzero: LaurentSeries,
    zero: LaurentSeries,
) -> LaurentSeries {
    let nonzero = nonzero.into_coefficients();
    let zero = zero.into_coefficients();
    LaurentSeries::new(
        if_nonzero_else(condition, nonzero[0].clone(), zero[0].clone()),
        if_nonzero_else(condition, nonzero[1].clone(), zero[1].clone()),
        if_nonzero_else(condition, nonzero[2].clone(), zero[2].clone()),
    )
}

pub(crate) fn finite_series(value: Atom) -> LaurentSeries {
    LaurentSeries::new(value, Atom::num(0), Atom::num(0))
}

pub(crate) fn triangle_one_mass(
    momenta: [&Atom; 3],
    mass: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let finite = finite_series(triangle_finite_one_mass(
        momenta[0], momenta[1], momenta[2], mass,
    ));
    let difference_2 = mass - momenta[1];
    let difference_3 = mass - momenta[2];
    let at_zero_first = if_series(
        &difference_2,
        if_series(
            &difference_3,
            triangle_ir_three(momenta[1], momenta[2], mass, mu_squared),
            triangle_ir_two(momenta[1], mass, mu_squared),
        ),
        if_series(
            &difference_3,
            triangle_ir_two(momenta[2], mass, mu_squared),
            triangle_ir_one(mass, mu_squared),
        ),
    );
    if_series(momenta[0], finite, at_zero_first)
}

pub(crate) fn triangle_two_masses(
    momenta: [&Atom; 3],
    mass_2: &Atom,
    mass_3: &Atom,
    mu_squared: &Atom,
) -> LaurentSeries {
    let finite = finite_series(triangle_finite_two_masses(
        momenta[0], momenta[1], momenta[2], mass_2, mass_3,
    ));
    let infrared = triangle_ir_four(momenta[1], mass_2, mass_3, mu_squared);
    let difference_1 = momenta[0] - mass_2;
    let difference_3 = momenta[2] - mass_3;
    if_series(
        &difference_1,
        finite.clone(),
        if_series(&difference_3, finite, infrared),
    )
}

pub(crate) fn with_triangle_normalization(series: LaurentSeries) -> LaurentSeries {
    let result = series.into_coefficients();
    let pi_squared = Symbol::PI.to_atom() * Symbol::PI.to_atom();
    LaurentSeries::new(
        &result[0] + &pi_squared * &result[2] / 12,
        result[1].clone(),
        result[2].clone(),
    )
}

/// Constructs a scalar triangle with its transparent native function definitions.
/// Momenta are [p1²,p2²,(p1+p2)²]; masses follow propagator order.
///
/// Exact quadratic/Gram-degenerate limits are not generally supported. In
/// particular, `p=[-4,-1,-1], m=[0,0,1]` currently returns an incorrect zero.
/// See STATUS.md for the development implementation's coverage limitations.
pub fn c0(p: [&Atom; 3], m: [&Atom; 3], mu_squared: &Atom) -> MappedLaurentSeries {
    let expressions = OneLoopExpressions::new();
    MappedLaurentSeries {
        series: expressions.c0(p, m, mu_squared),
        function_map: expressions.into_function_map(),
    }
}

/// Explicit alias for [c0].
pub fn c0_mapped(p: [&Atom; 3], m: [&Atom; 3], mu_squared: &Atom) -> MappedLaurentSeries {
    c0(p, m, mu_squared)
}

pub(crate) fn register_triangles(map: &mut FunctionMap) {
    let symbols =
        core::array::from_fn::<_, 7, _>(|i| symbolica::symbol!(format!("__olo_c0_arg_{i}")));
    let x = symbols.map(Symbol::to_atom);
    let p = [&x[0], &x[1], &x[2]];
    let m = [&x[3], &x[4], &x[5]];
    let lambda = (p[2] - p[0] - p[1]).pow(2) - Atom::num(4) * p[0] * p[1];
    let generic = if_nonzero_else(
        &sheet_exact::negative(&lambda),
        triangle_hv::expression(p, m),
        triangle_finite_three_masses(p, m),
    );
    let at_zero = triangle_vacuum(m);
    let nonzero_p = any_nonzero(p.into_iter().cloned());
    let definitions = [
        (0, triangle_massless(p, &x[6])),
        (1, triangle_one_mass(p, m[2], &x[6])),
        (3, triangle_two_masses(p, m[1], m[2], &x[6])),
        (
            7,
            finite_series(if_nonzero_else(&nonzero_p, generic, at_zero)),
        ),
    ];
    for (sector, series) in definitions {
        for (coefficient, body) in series.into_coefficients().into_iter().enumerate() {
            map.add_function_with_options(
                symbolica::symbol!(format!("__olo_c0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
                native_function_options(),
            )
            .expect("unique triangle sector");
        }
    }
    for (sector, canonical, shift) in [(2, 1, 2), (4, 1, 1), (5, 3, 1), (6, 3, 2)] {
        let args = (0..3)
            .map(|i| x[(i + shift) % 3].clone())
            .chain((0..3).map(|i| x[3 + (i + shift) % 3].clone()))
            .chain([x[6].clone()])
            .collect::<Vec<_>>();
        for coefficient in 0..3 {
            let body = symbolica::symbol!(format!("__olo_c0_sector_{canonical}_{coefficient}"))
                .call(&args);
            map.add_function(
                symbolica::symbol!(format!("__olo_c0_sector_{sector}_{coefficient}")),
                symbols.to_vec(),
                body,
            )
            .unwrap();
        }
    }
}

fn triangle_vacuum(m: [&Atom; 3]) -> Atom {
    let pair = |a: &Atom, b: &Atom, c: &Atom| {
        (log_over_one_minus(&(c / a)) - log_over_one_minus(&(c / b))) / (a - b)
    };
    let unequal = if_nonzero_else(
        &(m[0] - m[1]),
        pair(m[0], m[1], m[2]),
        pair(m[1], m[2], m[0]),
    );
    let difference = any_nonzero([m[0] - m[1], m[1] - m[2]]);
    if_nonzero_else(&difference, unequal, -Atom::num(1) / (Atom::num(2) * m[0]))
}
