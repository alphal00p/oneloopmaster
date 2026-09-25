//! Measure complete (unpruned) expression construction and subsequent selection.
use oneloop::{ExpressionOptions, ScalarIntegral, get_expression_with_options, select_branch};
use std::time::Instant;
use symbolica::prelude::*;

fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(run)
        .unwrap()
        .join()
        .unwrap();
}

fn run() {
    oneloop::initialize().unwrap();
    let args = std::env::args().collect::<Vec<_>>();
    let case = args.get(1).map(String::as_str).unwrap_or("C0-case");
    let budget = args
        .get(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(100_000_000);
    let a = if case == "C0-real" {
        symbol!("survey_real_a"; Real).to_atom()
    } else {
        symbol!("survey_a").to_atom()
    };
    let b = symbol!("survey_b").to_atom();
    let mu = symbol!("survey_mu"; Positive).to_atom();
    let (symbol, arguments, probes) = if case == "C0-case" || case == "C0-real" {
        (
            oneloop::C0(),
            vec![
                Atom::num(0),
                -&a,
                a.clone(),
                a.clone(),
                a.clone(),
                b.clone(),
                mu.clone(),
            ],
            vec![(a, 2), (b, 1), (mu, 1)],
        )
    } else {
        let (family, symbol) = match case {
            "A0" => (ScalarIntegral::A0, oneloop::A0()),
            "B0" => (ScalarIntegral::B0, oneloop::B0()),
            "dB0" => (ScalarIntegral::DB0, oneloop::dB0()),
            "C0" => (ScalarIntegral::C0, oneloop::C0()),
            "D0" => (ScalarIntegral::D0, oneloop::D0()),
            _ => panic!("unknown case"),
        };
        let arguments = (0..family.arity())
            .map(|i| symbol!(format!("survey_x{i}")).to_atom())
            .collect::<Vec<_>>();
        let momenta = match family {
            ScalarIntegral::A0 => 0,
            ScalarIntegral::B0 | ScalarIntegral::DB0 => 1,
            ScalarIntegral::C0 => 3,
            ScalarIntegral::D0 => 6,
        };
        let probes = arguments
            .iter()
            .enumerate()
            .map(|(i, x)| {
                (
                    x.clone(),
                    if i < momenta {
                        -(i as i32) - 1
                    } else {
                        i as i32 + 1
                    },
                )
            })
            .collect();
        (symbol, arguments, probes)
    };
    let start = Instant::now();
    eprintln!("starting {case}, max_nodes={budget}");
    if args.get(3).is_some_and(|s| s == "shared") {
        let series = oneloop::get_expression_shared_with_options(
            symbol.call(&arguments),
            ExpressionOptions {
                max_nodes: budget,
                max_depth: 512,
            },
        )
        .unwrap();
        eprintln!("shared construction: {:?}", start.elapsed());
        let rules = probes
            .into_iter()
            .map(|(x, n)| Replacement::new(x, n))
            .collect::<Vec<_>>();
        for (i, expression) in series.coefficients().iter().enumerate() {
            eprintln!(
                "coefficient {i}: {} bindings, {} bytes",
                expression.get_aliases().len(),
                expression.get_byte_size()
            );
            let start = Instant::now();
            let selected = oneloop::select_branch_shared(
                expression,
                &rules,
                ExpressionOptions {
                    max_nodes: budget,
                    max_depth: 4096,
                },
            )
            .unwrap();
            eprintln!(
                "late shared selection: {:?}, {} characters",
                start.elapsed(),
                selected.to_string().len()
            );
        }
        return;
    }
    let series = get_expression_with_options(
        symbol.call(&arguments),
        ExpressionOptions {
            max_nodes: budget,
            max_depth: 512,
        },
    );
    eprintln!("construction: {:?}", start.elapsed());
    let series = match series {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{e}");
            return;
        }
    };
    let rules = probes
        .into_iter()
        .map(|(x, n)| Replacement::new(x, n))
        .collect::<Vec<_>>();
    for (i, expression) in series.coefficients().iter().enumerate() {
        let mut nodes = 0usize;
        let mut branches = 0usize;
        expression.visitor(&mut |v| {
            nodes += 1;
            branches += usize::from(v.get_symbol() == Some(Symbol::IF));
            true
        });
        eprintln!(
            "coefficient {i}: {nodes} nodes, {branches} ifs, {} bytes",
            expression.as_view().get_byte_size()
        );
        let start = Instant::now();
        let selected = select_branch(expression, &rules);
        eprintln!(
            "selection: {:?}, {} characters",
            start.elapsed(),
            selected.to_string().len()
        );
    }
}
