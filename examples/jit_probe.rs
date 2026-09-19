//! Development probe for native SymJIT compilation; not a performance claim.
use std::time::Instant;
use symbolica::prelude::*;

fn main() {
    let family = std::env::args().nth(1).unwrap_or_else(|| "A0".into());
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || run(&family))
        .unwrap()
        .join()
        .unwrap();
}

fn run(family: &str) {
    let arity = match family {
        "A0" => 2,
        "B0" | "dB0" => 4,
        "C0" => 7,
        "D0" => 11,
        _ => panic!("expected A0, B0, dB0, C0, or D0"),
    };
    let args: Vec<_> = (0..arity)
        .map(|i| symbol!(format!("jit_probe_{i}")).to_atom())
        .collect();
    let context = oneloop::OneLoopExpressions::new();
    let series = match family {
        "A0" => oneloop::a0(&args[0], &args[1]),
        "B0" => oneloop::b0(&args[0], &args[1], &args[2], &args[3]),
        "dB0" => oneloop::db0(&args[0], &args[1], &args[2], &args[3]),
        "C0" => context.c0(
            [&args[0], &args[1], &args[2]],
            [&args[3], &args[4], &args[5]],
            &args[6],
        ),
        "D0" => context.d0(
            core::array::from_fn(|i| &args[i]),
            core::array::from_fn(|i| &args[i + 6]),
            &args[10],
        ),
        _ => unreachable!(),
    };
    let start = Instant::now();
    let exact = context.evaluator(series.coefficients(), &args).unwrap();
    eprintln!("{family} native build: {:?}", start.elapsed());
    let start = Instant::now();
    let mut jit = context.jit_evaluator(series.coefficients(), &args).unwrap();
    eprintln!(
        "{family} JIT O2: {:?}; cache bytes: {}",
        start.elapsed(),
        jit.to_bytes().unwrap().len()
    );
    let input = match family {
        "A0" => vec![2., 1.],
        "B0" | "dB0" => vec![-1., 1., 2., 1.],
        "C0" => vec![-1., -2., -3., 1., 2., 3., 1.],
        "D0" => vec![-1., -2., -3., -4., -12., -15., 1., 2., 3., 4., 1.],
        _ => unreachable!(),
    }
    .into_iter()
    .map(|v| Complex::new(v, 0.))
    .collect::<Vec<_>>();
    let mut native = exact.map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let mut expected = [Complex::new(0., 0.); 3];
    let mut actual = expected;
    native.evaluate(&input, &mut expected);
    jit.evaluate(&input, &mut actual).unwrap();
    eprintln!("native {expected:?}; JIT {actual:?}");
    for (a, b) in actual.into_iter().zip(expected) {
        let error = (a.re - b.re).hypot(a.im - b.im);
        assert!(error.is_finite() && error < 1e-10);
    }
}
