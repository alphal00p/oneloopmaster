//! Warm scalar and row-major batch throughput against shared Fortran fixtures.
//! Usage: performance_survey [calls=20000] [repetitions=7] [family=all] [route=jit] [batches=1,4,32,256,1024]
//! `native` uses ahead-of-time generic Rust; `jit` uses the prepared SymJIT backend;
//! `symbols` compiles public master
//! Symbols with their transparent native FunctionMap and supports the same batches.
//! `hooks-jit` compiles the bare public calls without a FunctionMap and supports
//! the same batches, invoking the automatically registered numerical callbacks.
//! `hooks` keeps that bare outer evaluator interpreted and is scalar-only.
//! Compare interpreted `hooks` TSV only with --batch 1.
//! Emits TSV; construction, correctness checks, and allocation are outside timing.
#[path = "../tests/support/fixtures.rs"]
mod fixtures;
use fixtures::Fixture;
use oneloop::{JitEvaluator, NativeEvaluator, OneLoopExpressions, ScalarEvaluator, ScalarIntegral};
use std::{hint::black_box, time::Instant};
use symbolica::{
    evaluate::{ExpressionEvaluator, JITCompiledEvaluator},
    prelude::*,
};

#[derive(Clone, Copy, Debug)]
enum Route {
    Native,
    Jit,
    Symbols,
    Hooks,
    HooksJit,
}

impl SurveyEvaluator for NativeEvaluator<f64> {
    fn family(&self) -> ScalarIntegral {
        NativeEvaluator::family(self)
    }
    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) {
        self.evaluate(input, output).unwrap();
    }
    fn batch(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>], rows: usize) {
        self.evaluate_batch(input, output, rows).unwrap();
    }
}

trait SurveyEvaluator {
    fn family(&self) -> ScalarIntegral;
    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]);
    fn batch(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>], rows: usize);
}

impl SurveyEvaluator for ScalarEvaluator {
    fn family(&self) -> ScalarIntegral {
        ScalarEvaluator::family(self)
    }

    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) {
        self.evaluate(input, output).unwrap();
    }

    fn batch(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>], rows: usize) {
        self.evaluate_batch(input, output, rows).unwrap();
    }
}

fn master_expressions(family: ScalarIntegral) -> ([Atom; 3], Vec<Atom>) {
    let master = match family {
        ScalarIntegral::A0 => oneloop::A0(),
        ScalarIntegral::B0 => oneloop::B0(),
        ScalarIntegral::DB0 => oneloop::dB0(),
        ScalarIntegral::C0 => oneloop::C0(),
        ScalarIntegral::D0 => oneloop::D0(),
    };
    let parameters = (0..family.arity())
        .map(|index| symbol!(format!("survey_hook_{}_{}", family.name(), index)).to_atom())
        .collect::<Vec<_>>();
    let expressions = [0, -1, -2].map(|tag| {
        let mut arguments = vec![Atom::num(tag)];
        arguments.extend_from_slice(&parameters);
        master.call(&arguments)
    });
    (expressions, parameters)
}

struct SymbolEvaluator {
    family: ScalarIntegral,
    inner: JitEvaluator,
}

impl SymbolEvaluator {
    fn new(family: ScalarIntegral) -> Self {
        let (expressions, parameters) = master_expressions(family);
        // The public calls remain the input expressions, but native definitions
        // take precedence over their automatic scalar callbacks during compilation.
        let inner = OneLoopExpressions::new()
            .jit_evaluator(&expressions, &parameters)
            .expect("compile public master Symbols with the native FunctionMap");
        Self { family, inner }
    }
}

impl SurveyEvaluator for SymbolEvaluator {
    fn family(&self) -> ScalarIntegral {
        self.family
    }

    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) {
        self.inner.evaluate(input, output).unwrap();
    }

    fn batch(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>], rows: usize) {
        self.inner.evaluate_batch(input, output, rows).unwrap();
    }
}

struct HookEvaluator {
    family: ScalarIntegral,
    inner: ExpressionEvaluator<Complex<f64>>,
}

impl HookEvaluator {
    fn new(family: ScalarIntegral) -> Self {
        let (expressions, parameters) = master_expressions(family);
        // Deliberately omit the native FunctionMap: this measures the public
        // scalar callback path and its cache locking, not expanded definitions.
        let inner = Atom::evaluator_multiple(&expressions, &parameters)
            .direct_translation(true)
            .build()
            .unwrap()
            .map_coeff(&|value| Complex::new(value.re.to_f64(), value.im.to_f64()));
        Self { family, inner }
    }
}

impl SurveyEvaluator for HookEvaluator {
    fn family(&self) -> ScalarIntegral {
        self.family
    }

    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) {
        self.inner.evaluate(input, output);
    }

    fn batch(&mut self, _input: &[Complex<f64>], _output: &mut [Complex<f64>], _rows: usize) {
        unreachable!("public Symbol hooks are benchmarked only at batch size 1");
    }
}

struct HookJitEvaluator {
    family: ScalarIntegral,
    inner: JITCompiledEvaluator<Complex<f64>>,
}

impl HookJitEvaluator {
    fn new(family: ScalarIntegral) -> Self {
        let (expressions, parameters) = master_expressions(family);
        // No native FunctionMap: the outer JIT retains the public master
        // callbacks, whose automatically prepared backends do the numerical work.
        let exact = Atom::evaluator_multiple(&expressions, &parameters)
            .direct_translation(true)
            .build()
            .expect("build bare public master expressions");
        let inner = exact
            .jit_compile::<Complex<f64>>(oneloop::jit_settings())
            .expect("JIT-compile bare public master callbacks");
        assert_eq!(inner.input_count(), family.arity());
        assert_eq!(inner.output_count(), 3);
        Self { family, inner }
    }
}

impl SurveyEvaluator for HookJitEvaluator {
    fn family(&self) -> ScalarIntegral {
        self.family
    }

    fn point(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>]) {
        self.inner.evaluate(input, output);
    }

    fn batch(&mut self, input: &[Complex<f64>], output: &mut [Complex<f64>], rows: usize) {
        self.inner.batch_evaluate(input, output, rows);
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let calls: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(20000);
    let repetitions: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(7);
    let family = args.next().filter(|s| !s.eq_ignore_ascii_case("all"));
    let route = match args.next().as_deref().unwrap_or("jit") {
        value if value.eq_ignore_ascii_case("native") => Route::Native,
        value if value.eq_ignore_ascii_case("jit") => Route::Jit,
        value if value.eq_ignore_ascii_case("symbols") => Route::Symbols,
        value if value.eq_ignore_ascii_case("hooks") => Route::Hooks,
        value if value.eq_ignore_ascii_case("hooks-jit") => Route::HooksJit,
        value => panic!("unknown route {value}; expected native, jit, symbols, hooks or hooks-jit"),
    };
    let batches = args.next().map(|value| {
        let batches = value
            .split(',')
            .map(|size| size.parse::<usize>().expect("positive batch size"))
            .collect::<Vec<_>>();
        assert!(!batches.is_empty() && batches.iter().all(|&size| size > 0));
        let mut unique = batches.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), batches.len(), "duplicate batch size");
        if matches!(route, Route::Hooks) {
            assert_eq!(batches, [1], "public Symbol hooks are scalar-only");
        }
        batches
    });
    assert!(
        family
            .as_ref()
            .is_none_or(|s| ["A0", "B0", "dB0", "C0", "D0"]
                .iter()
                .any(|name| s.eq_ignore_ascii_case(name))),
        "unknown scalar family"
    );
    assert!(calls > 0 && repetitions > 0 && args.next().is_none());
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || run(calls, repetitions, family, route, batches))
        .unwrap()
        .join()
        .unwrap();
}

fn run(
    calls: usize,
    repetitions: usize,
    selected: Option<String>,
    route: Route,
    batches: Option<Vec<usize>>,
) {
    let startup = Instant::now();
    oneloop::initialize().unwrap();
    eprintln!("all-family eager initialization {:?}", startup.elapsed());
    let fixtures = fixtures::parse(include_str!("../tests/data/benchmark.txt"));
    eprintln!(
        "survey route={route:?}; all three Laurent outputs; warmup=256; interpreted hooks are scalar-only"
    );
    println!("mode\tfamily\tname\tbatch\tcalls\trepetition\tns_per_call\tchecksum_re\tchecksum_im");
    for family in [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ] {
        if selected
            .as_ref()
            .is_some_and(|s| !s.eq_ignore_ascii_case(family.name()))
        {
            continue;
        }
        let started = Instant::now();
        let rows: Vec<_> = fixtures.iter().filter(|r| r.family == family).collect();
        match route {
            Route::Native => {
                let evaluator = NativeEvaluator::<f64>::new(family).unwrap();
                eprintln!(
                    "route=native {} exact-constant setup {:?}",
                    family.name(),
                    started.elapsed()
                );
                family_workloads(
                    evaluator,
                    &rows,
                    calls,
                    repetitions,
                    batches.as_deref().unwrap_or(&[1, 4, 32, 256, 1024]),
                );
            }
            Route::Jit => {
                #[cfg(feature = "prebuilt")]
                let evaluator = ScalarEvaluator::prebuilt(family).unwrap();
                #[cfg(not(feature = "prebuilt"))]
                let evaluator = ScalarEvaluator::rebuild(family).unwrap();
                eprintln!(
                    "route=jit {} prepared-backend clone {:?}",
                    family.name(),
                    started.elapsed()
                );
                family_workloads(
                    evaluator,
                    &rows,
                    calls,
                    repetitions,
                    batches.as_deref().unwrap_or(&[1, 4, 32, 256, 1024]),
                );
            }
            Route::Symbols => {
                let evaluator = SymbolEvaluator::new(family);
                eprintln!(
                    "route=symbols {} public master expressions + native FunctionMap JIT compilation {:?}",
                    family.name(),
                    started.elapsed()
                );
                family_workloads(
                    evaluator,
                    &rows,
                    calls,
                    repetitions,
                    batches.as_deref().unwrap_or(&[1, 4, 32, 256, 1024]),
                );
            }
            Route::Hooks => {
                let evaluator = HookEvaluator::new(family);
                eprintln!(
                    "route=hooks {} expression build {:?}; using already prepared shared backend",
                    family.name(),
                    started.elapsed()
                );
                family_workloads(evaluator, &rows, calls, repetitions, &[1]);
            }
            Route::HooksJit => {
                let evaluator = HookJitEvaluator::new(family);
                eprintln!(
                    "route=hooks-jit {} bare public master JIT compilation {:?}; using automatic cached callbacks without a FunctionMap",
                    family.name(),
                    started.elapsed()
                );
                family_workloads(
                    evaluator,
                    &rows,
                    calls,
                    repetitions,
                    batches.as_deref().unwrap_or(&[1, 4, 32, 256, 1024]),
                );
            }
        }
    }
}

fn family_workloads<E: SurveyEvaluator>(
    mut evaluator: E,
    rows: &[&Fixture],
    calls: usize,
    repetitions: usize,
    batches: &[usize],
) {
    for row in rows {
        let mut output = [Complex::new(f64::NAN, f64::NAN); 3];
        evaluator.point(&row.args, &mut output);
        row.check(&output);
    }
    for &batch in batches {
        for row in rows {
            workload(
                &mut evaluator,
                &[*row],
                "SAME",
                &row.name,
                calls,
                repetitions,
                batch,
            );
        }
        workload(
            &mut evaluator,
            rows,
            "HETERO",
            "all",
            calls,
            repetitions,
            batch,
        );
    }
}

fn workload<E: SurveyEvaluator>(
    evaluator: &mut E,
    rows: &[&Fixture],
    mode: &str,
    name: &str,
    calls: usize,
    repetitions: usize,
    batch: usize,
) {
    let arity = evaluator.family().arity();
    let inputs: Vec<_> = rows
        .iter()
        .cycle()
        .take(calls)
        .flat_map(|r| r.args.iter().copied())
        .collect();
    let mut outputs = vec![Complex::new(f64::NAN, f64::NAN); calls * 3];
    // Verify the precise mixed-lane layout, including a partial last batch.
    evaluate_all(evaluator, &inputs, &mut outputs, arity, batch);
    for (index, output) in outputs.chunks_exact(3).enumerate() {
        rows[index % rows.len()].check(output);
    }
    let mut warm = [Complex::new(0., 0.); 3];
    for row in rows.iter().cycle().take(256) {
        evaluator.point(&row.args, &mut warm);
    }
    for repetition in 1..=repetitions {
        let started = Instant::now();
        let checksum = evaluate_all(evaluator, black_box(&inputs), &mut outputs, arity, batch);
        let elapsed = started.elapsed().as_nanos() as f64 / calls as f64;
        black_box(&outputs);
        println!(
            "{mode}\t{}\t{name}\t{batch}\t{calls}\t{repetition}\t{elapsed:.9}\t{:.17e}\t{:.17e}",
            evaluator.family().name(),
            checksum.re,
            checksum.im
        );
    }
}

fn evaluate_all<E: SurveyEvaluator>(
    evaluator: &mut E,
    inputs: &[Complex<f64>],
    outputs: &mut [Complex<f64>],
    arity: usize,
    batch: usize,
) -> Complex<f64> {
    let mut checksum = Complex::new(0., 0.);
    for (input, output) in inputs
        .chunks(arity * batch)
        .zip(outputs.chunks_mut(3 * batch))
    {
        if batch == 1 {
            evaluator.point(input, output);
        } else {
            evaluator.batch(input, output, input.len() / arity);
        }
        for row in output.chunks_exact(3) {
            checksum += row[0];
        }
    }
    checksum
}
