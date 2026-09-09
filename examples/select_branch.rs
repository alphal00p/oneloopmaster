//! Print the full B0 expression and a parametric formula for one analytic region.
use oneloop::{B0, get_expression, select_branch};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::Float,
    id::Replacement,
    printer::PrintOptions,
};

fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let psq = symbolica::symbol!("oneloop_example::psq"; Real).to_atom();
            let m2 = symbolica::symbol!("oneloop_example::m2"; Positive).to_atom();
            let master = B0().call((&psq, &m2, &m2, 1));
            let series = get_expression(master).expect("complete B0 expression");
            let rules = [
                Replacement::new(psq, Atom::num(Float::parse("3.23", Some(128)).unwrap())),
                Replacement::new(m2, 1),
            ];
            for (tag, coefficient) in [0, -1, -2].into_iter().zip(series.coefficients()) {
                println!("# all_branches, epsilon^{tag}");
                println!("{}", coefficient.printer(PrintOptions::full()));
                println!("# selected at psq=3.23, m2=1, epsilon^{tag}");
                println!(
                    "{}",
                    select_branch(coefficient, &rules).printer(PrintOptions::full())
                );
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
