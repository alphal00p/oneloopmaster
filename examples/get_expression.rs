//! Print complete A0 and massless C0 bodies. No numerical evaluator is built here.
use oneloop::{ScalarIntegral, get_expression};
use symbolica::{
    atom::{Atom, AtomCore},
    printer::PrintOptions,
};

fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let mass = symbolica::symbol!("example_mass_squared"; Positive).to_atom();
            let scale = symbolica::symbol!("example_mu_squared"; Positive).to_atom();
            let invariant = symbolica::symbol!("example_s_squared"; Real).to_atom();
            for (family, master, arguments) in [
                (ScalarIntegral::A0, oneloop::A0(), vec![mass, scale.clone()]),
                (
                    ScalarIntegral::C0,
                    oneloop::C0(),
                    vec![
                        Atom::num(0),
                        Atom::num(0),
                        invariant,
                        Atom::num(0),
                        Atom::num(0),
                        Atom::num(0),
                        scale,
                    ],
                ),
            ] {
                let series = get_expression(master.call(&arguments))
                    .expect("complete expression within budget");
                for (tag, body) in [0, -1, -2].into_iter().zip(series.coefficients()) {
                    println!(
                        "{}[{tag}] = {}",
                        family.name(),
                        body.printer(PrintOptions::full())
                    );
                }
            }
        })
        .unwrap()
        .join()
        .unwrap();
}
