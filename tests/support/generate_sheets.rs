//! Optional fixture generator; compile with rustc, not cargo test.
//! Requires the sibling standalone numerical source checkout. Normal tests do not.
#[path = "../../../symbolica-oneloop/src/complex.rs"]
mod complex;
pub use complex::Complex64;
#[path = "../../../symbolica-oneloop/src/dilog.rs"]
mod dilog;
#[path = "../../../symbolica-oneloop/src/math.rs"]
mod math;
fn main() {
    let points = [
        (0.2, 0.),
        (2., 0.),
        (-2., 0.),
        (-0.2, 0.),
        (0.2, 0.3),
        (-0.2, 0.3),
        (-0.2, -0.3),
        (2., 3.),
        (-2., -3.),
    ];
    for (a, b) in points {
        for (c, d) in points {
            for sign in [-1, 1] {
                for div in [false, true] {
                    let x = math::SheetComplex::from_value(Complex64::new(a, b), sign);
                    let y = math::SheetComplex::from_value(Complex64::new(c, d), -sign);
                    let z = if div { x / y } else { x * y };
                    let log = z.log();
                    let li = dilog::dilog(z);
                    println!(
                        "{a} {b} {c} {d} {sign} {} {:.17e} {:.17e} {} {:.17e} {:.17e} {:.17e} {:.17e}",
                        u8::from(div),
                        z.value().re,
                        z.value().im,
                        z.phase(),
                        log.re,
                        log.im,
                        li.re,
                        li.im
                    );
                }
            }
        }
    }
}
