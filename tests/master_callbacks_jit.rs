//! Bare public-master outer JIT: automatic callbacks, never a native FunctionMap.
#[path = "support/fixtures.rs"]
mod fixtures;

use fixtures::Fixture;
use oneloop::{ScalarEvaluator, ScalarIntegral};
use symbolica::{evaluate::JITCompiledEvaluator, prelude::*};

type C = Complex<f64>;

fn compile_bare(family: ScalarIntegral, order: [usize; 3]) -> JITCompiledEvaluator<C> {
    let master = match family {
        ScalarIntegral::A0 => oneloop::A0(),
        ScalarIntegral::B0 => oneloop::B0(),
        ScalarIntegral::DB0 => oneloop::dB0(),
        ScalarIntegral::C0 => oneloop::C0(),
        ScalarIntegral::D0 => oneloop::D0(),
    };
    let parameters = (0..family.arity())
        .map(|index| symbol!(format!("bare_master_jit_{}_{}", family.name(), index)).to_atom())
        .collect::<Vec<_>>();
    let expressions = order.map(|coefficient| {
        let mut arguments = vec![Atom::num([0, -1, -2][coefficient])];
        // Numeric arguments are momenta, squared masses, then mu_squared.
        arguments.extend_from_slice(&parameters);
        master.call(&arguments)
    });
    // Deliberately do not attach OneLoopExpressions or any FunctionMap. The
    // compiled outer graph must retain the automatic public-master callbacks.
    let exact = Atom::evaluator_multiple(&expressions, &parameters)
        .direct_translation(true)
        .build()
        .expect("build bare public-master expressions");
    let evaluator = exact
        .jit_compile::<C>(oneloop::jit_settings())
        .expect("compile bare public-master callbacks");
    assert!(evaluator.has_external_functions());
    assert_eq!(evaluator.input_count(), family.arity());
    assert_eq!(evaluator.output_count(), 3);
    evaluator
}

fn record_comparison(
    actual: C,
    expected: C,
    normalization: f64,
    label: &str,
    failures: &mut Vec<String>,
) {
    // Same dimension normalization and tolerances as tests/masters.rs.
    let error = (actual.re - expected.re).hypot(actual.im - expected.im) * normalization;
    let tolerance = 2e-10 + 2e-8 * expected.re.hypot(expected.im) * normalization;
    if !actual.re.is_finite() || !actual.im.is_finite() || !error.is_finite() || error > tolerance {
        failures.push(format!(
            "{label}: actual={actual:?}, expected={expected:?}, normalized error={error}, tolerance={tolerance}"
        ));
    }
}

fn check_row(
    row: &Fixture,
    prepared: &[C; 3],
    actual: &[C],
    order: [usize; 3],
    label: &str,
    failures: &mut Vec<String>,
) {
    assert_eq!(actual.len(), 3);
    let power = match row.family {
        ScalarIntegral::A0 => -1,
        ScalarIntegral::B0 => 0,
        ScalarIntegral::DB0 | ScalarIntegral::C0 => 1,
        ScalarIntegral::D0 => 2,
    };
    let normalization = row.args.last().unwrap().re.powi(power);
    for (position, coefficient) in order.into_iter().enumerate() {
        let prefix = format!(
            "{} {} {label}, tag {} (output {position})",
            row.family.name(),
            row.name,
            [0, -1, -2][coefficient]
        );
        record_comparison(
            actual[position],
            row.expected[coefficient],
            normalization,
            &format!("{prefix} versus fixture"),
            failures,
        );
        record_comparison(
            actual[position],
            prepared[coefficient],
            normalization,
            &format!("{prefix} versus prepared backend"),
            failures,
        );
    }
}

fn exercise(
    evaluator: &mut JITCompiledEvaluator<C>,
    rows: &[&Fixture],
    prepared: &[[C; 3]],
    indices: &[usize],
    order: [usize; 3],
    stage: &str,
    failures: &mut Vec<String>,
) {
    for &index in indices {
        let mut output = [C::new(f64::NAN, f64::NAN); 3];
        evaluator.evaluate(&rows[index].args, &mut output);
        check_row(
            rows[index],
            &prepared[index],
            &output,
            order,
            &format!("{stage} scalar, order={order:?}"),
            failures,
        );
    }
    for batch in [1, 3, 4, 5, 31] {
        for chunk in indices.chunks(batch) {
            let input: Vec<_> = chunk
                .iter()
                .flat_map(|&index| rows[index].args.iter().copied())
                .collect();
            let mut output = vec![C::new(f64::NAN, f64::NAN); 3 * chunk.len()];
            evaluator.batch_evaluate(&input, &mut output, chunk.len());
            for (&index, actual) in chunk.iter().zip(output.chunks_exact(3)) {
                check_row(
                    rows[index],
                    &prepared[index],
                    actual,
                    order,
                    &format!("{stage} batch={batch}, order={order:?}"),
                    failures,
                );
            }
        }
    }
}

#[test]
fn bare_master_outer_jit_and_portable_batches_match_all_acceptance_rows() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            oneloop::initialize().unwrap();
            let fixtures = fixtures::parse(include_str!("data/parity.txt"));
            let families = [
                ScalarIntegral::A0,
                ScalarIntegral::B0,
                ScalarIntegral::DB0,
                ScalarIntegral::C0,
                ScalarIntegral::D0,
            ];
            let orders = [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            let mut inventory = [0; 5];
            let mut failures = Vec::new();
            for family in families {
                let rows: Vec<_> = fixtures.iter().filter(|row| row.family == family).collect();
                inventory[family as usize] = rows.len();
                assert!(rows.iter().all(|row| row.args.len() == family.arity()));
                let mut prepared = ScalarEvaluator::cached(family).unwrap();
                let reference: Vec<_> = rows
                    .iter()
                    .map(|row| {
                        let mut output = [C::new(f64::NAN, f64::NAN); 3];
                        prepared.evaluate(&row.args, &mut output).unwrap();
                        output
                    })
                    .collect();
                let pole = rows
                    .iter()
                    .position(|row| row.expected[1..].iter().any(|value| *value != C::new(0., 0.)))
                    .expect("every acceptance family includes an infrared or ultraviolet pole");
                // Every acceptance row appears twice adjacently. This exercises
                // identical SIMD lanes as well as transitions to other points.
                // The final three rows add a nonadjacent return and partial tails.
                let all_rows: Vec<_> = (0..rows.len())
                    .flat_map(|index| [index, index])
                    .chain([0, rows.len() - 1, 0])
                    .collect();
                let controls = [0, pole, rows.len() - 1, 0, pole, rows.len() - 1, 0, 0];
                for order in orders {
                    let mut original = compile_bare(family, order);
                    let bytes = original.export_portable().unwrap();
                    let mut restored =
                        JITCompiledEvaluator::<C>::import_portable(&bytes, oneloop::jit_settings())
                            .expect("restore bare public-master callback definitions");
                    assert!(restored.has_external_functions());
                    assert_eq!(restored.input_count(), family.arity());
                    assert_eq!(restored.output_count(), 3);
                    // The canonical order covers all 293 rows; each other tag
                    // permutation covers repeated/mixed finite and pole controls.
                    let indices = if order == [0, 1, 2] {
                        all_rows.as_slice()
                    } else {
                        controls.as_slice()
                    };
                    for (stage, evaluator) in
                        [("original", &mut original), ("restored", &mut restored)]
                    {
                        exercise(
                            evaluator,
                            &rows,
                            &reference,
                            indices,
                            order,
                            stage,
                            &mut failures,
                        );
                    }
                }
                eprintln!(
                    "{}: {} acceptance rows, scalar + five batch layouts, six tag orders, original/restored bare JIT",
                    family.name(),
                    rows.len()
                );
            }
            assert_eq!(inventory, [4, 13, 12, 125, 139]);
            assert_eq!(inventory.iter().sum::<usize>(), 293);
            assert!(
                failures.is_empty(),
                "{} bare-master outer-JIT coefficient comparisons failed:\n{}",
                failures.len(),
                failures.join("\n")
            );
        })
        .unwrap()
        .join()
        .unwrap();
}
