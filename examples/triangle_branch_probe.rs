//! Regression diagnostic of native/JIT one-mass triangle sheet intermediates.
use symbolica::evaluate::ExpressionEvaluator;
use symbolica::prelude::*;

fn call(name: &str, args: &[Atom]) -> Atom {
    symbol!(format!("oneloop::{name}")).call(args)
}
fn choose(c: &Atom, yes: Atom, no: Atom) -> Atom {
    Symbol::IF.call((c, yes, no))
}
fn real(z: &Atom) -> Atom {
    call("__olo_real_part", std::slice::from_ref(z))
}
fn imag(z: &Atom) -> Atom {
    call("__olo_imaginary_part", std::slice::from_ref(z))
}
fn sign(z: &Atom) -> Atom {
    call("__olo_sign_nonnegative", std::slice::from_ref(z))
}
fn negative(z: &Atom) -> Atom {
    (1 - sign(z)) / 2
}
fn sqrt(z: &Atom) -> Atom {
    call("__olo_sqrt_lower", std::slice::from_ref(z))
}
#[derive(Clone)]
struct Q {
    v: Atom,
    p: Atom,
    o: Atom,
}
impl Q {
    fn new(v: Atom, upper: bool) -> Self {
        let n = negative(&real(&v));
        let im = imag(&v);
        let s = choose(&im, sign(&im), Atom::num(if upper { 1 } else { -1 }));
        Self {
            v: (1 - 2 * &n) * v,
            p: &n * s,
            o: n,
        }
    }
    fn product(&self, b: &Self, divide: bool) -> Self {
        let args = [
            self.v.clone(),
            self.p.clone(),
            self.o.clone(),
            b.v.clone(),
            b.p.clone(),
            b.o.clone(),
        ];
        let name = if divide { "quotient" } else { "product" };
        Self {
            v: call(&format!("__olo_sheet_{name}_0"), &args),
            p: call(&format!("__olo_sheet_{name}_1"), &args),
            o: call(&format!("__olo_sheet_{name}_2"), &args),
        }
    }
    fn scale(&self, s: &Atom) -> Self {
        Self {
            v: &self.v * s,
            p: self.p.clone(),
            o: self.o.clone(),
        }
    }
    fn value(&self) -> Atom {
        (1 - 2 * &self.o) * &self.v
    }
    fn log(&self) -> Atom {
        call(
            "__olo_sheet_log",
            &[self.v.clone(), self.p.clone(), self.o.clone()],
        )
    }
    fn dilog(&self) -> Atom {
        call(
            "__olo_sheet_dilog",
            &[self.v.clone(), self.p.clone(), self.o.clone()],
        )
    }
    fn logc(&self) -> Atom {
        let d = 1 - self.value();
        choose(&d, self.log() / &d, Atom::num(-1))
    }
}
fn dd(a: &Q, b: &Q) -> Atom {
    let d = a.value() - b.value();
    choose(&d, (a.dilog() - b.dilog()) / &d, a.logc())
}
fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(run)
        .unwrap()
        .join()
        .unwrap();
}
fn run() {
    let x = (0..7)
        .map(|i| symbol!(format!("triangle_probe_{i}")).to_atom())
        .collect::<Vec<_>>();
    let (p1, p2, p3, m) = (&x[0], &x[1], &x[2], &x[5]);
    let sm = sqrt(m);
    let sr = (m * m.conj()).sqrt().sqrt();
    let r23 = -p1 / (&sr * &sr);
    let r24 = (m - p3) / (&sr * &sm);
    let r34 = (m - p2) / (&sr * &sm);
    let a = &r34 * &r24 - &r23;
    let b = &r24 / &sr + &r34 / &sr - &r23 / &sm;
    let c = Atom::num(1) / (&sr * &sr);
    let d = &b * &b - 4 * &a * &c;
    let sd = sqrt(&d);
    let r1: Atom = (&b + &sd) / (2 * &a);
    let r2: Atom = (&b - &sd) / (2 * &a);
    let q1 = Q::new(r1.clone(), true);
    let q2 = Q::new(r2.clone(), true);
    let qm = Q::new(sm.clone(), false);
    let q23 = Q::new(r23.clone(), false);
    let q24 = Q::new(r24.clone(), false);
    let q34 = Q::new(r34.clone(), false);
    let ratio = q1.product(&q2, true);
    let rp = q1.product(&q2, false);
    let mp = qm.product(&qm, false);
    let rm = rp.product(&mp, true);
    let q1m = q1.product(&qm, false);
    let q2m = q2.product(&qm, false);
    let q1a = q1.product(&q24.scale(&sr), false);
    let q2a = q2.product(&q24.scale(&sr), false);
    let q1b = q1.product(&q34.scale(&sr), false);
    let q2b = q2.product(&q34.scale(&sr), false);
    let context = oneloop::OneLoopExpressions::new();
    let full = context.c0([p1, p2, p3], [&x[3], &x[4], m], &x[6]);
    let term_log = -(rm.log() / 2 + q23.scale(&(&sr * &sr)).log()) * ratio.logc() / &r2;
    let term_m = -dd(&q1m, &q2m) * &sm;
    let term24 = dd(&q1a, &q2a) * &r24 * &sr;
    let term34 = dd(&q1b, &q2b) * &r34 * &sr;
    let finite = (&term_log + &term_m + &term24 + &term34) / (&a * &sr * &sr * &sm);
    let mut labels = vec!["full".to_string(), "manual".into()];
    let mut parts = vec![full.coefficients()[0].clone(), finite];
    for (label, value) in [
        ("norm_square", m * m.conj()),
        ("sqrt_mass", sm),
        ("sqrt_real", sr),
        ("r23", r23),
        ("r24", r24),
        ("r34", r34),
        ("a", a),
        ("b", b),
        ("c", c),
        ("d", d),
        ("sqrt_d", sd),
        ("root1", r1),
        ("root2", r2),
        ("term_log", term_log),
        ("term_m", term_m),
        ("term24", term24),
        ("term34", term34),
    ] {
        labels.push(label.into());
        parts.push(value);
    }
    for (name, q) in [
        ("q1", q1),
        ("q2", q2),
        ("qm", qm),
        ("ratio", ratio),
        ("rp", rp),
        ("mp", mp),
        ("rm", rm),
        ("q1m", q1m),
        ("q2m", q2m),
        ("q1a", q1a),
        ("q2a", q2a),
        ("q1b", q1b),
        ("q2b", q2b),
    ] {
        for (component, value) in [
            ("v", q.v.clone()),
            ("p", q.p.clone()),
            ("o", q.o.clone()),
            ("log", q.log()),
            ("li2", q.dilog()),
        ] {
            labels.push(format!("{name}.{component}"));
            parts.push(value);
        }
    }
    let exact = context.evaluator(&parts, &x).unwrap();
    let mut native = exact
        .clone()
        .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
    let mut jit = native.jit_compile(oneloop::jit_settings()).unwrap();
    let bytes = bincode::encode_to_vec(&native, bincode::config::standard()).unwrap();
    let (source, consumed): (ExpressionEvaluator<Complex<f64>>, usize) =
        bincode::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
    assert_eq!(consumed, bytes.len());
    let mut restored = source.jit_compile(oneloop::jit_settings()).unwrap();
    let mut failures = 0;
    for (line, values) in [
        (
            100,
            [3., 1., 2., 1.581_143_556_845_239_2, -0.07289733030227923],
        ),
        (107, [2., 3., -1., 0.5656924248442534, -0.08859878061223221]),
    ] {
        let input = [
            Complex::new(values[0], 0.),
            Complex::new(values[1], 0.),
            Complex::new(values[2], 0.),
            Complex::new(0., 0.),
            Complex::new(0., 0.),
            Complex::new(values[3], values[4]),
            Complex::new(1., 0.),
        ];
        let mut expected = vec![Complex::new(f64::NAN, f64::NAN); parts.len()];
        native.evaluate(&input, &mut expected);
        for (stage, compiled) in [("original", &mut jit), ("restored", &mut restored)] {
            let mut actual = vec![Complex::new(f64::NAN, f64::NAN); parts.len()];
            compiled.evaluate(&input, &mut actual);
            for ((label, native), jit) in labels.iter().zip(&expected).zip(&actual) {
                let delta = (native.re - jit.re).hypot(native.im - jit.im);
                if !delta.is_finite()
                    || delta > 1e-11 * (1. + native.re.hypot(native.im))
                    || (label == "norm_square" && jit.im != 0.)
                {
                    failures += 1;
                    eprintln!("line{line}/{stage} {label}: native={native:?} JIT={jit:?}");
                }
            }
            println!(
                "line{line}/{stage}: native={:?}, JIT={:?}, norm_square={:?}",
                expected[0], actual[0], actual[2]
            );
        }
    }
    assert_eq!(failures, 0, "native/JIT triangle intermediate mismatches");
}
