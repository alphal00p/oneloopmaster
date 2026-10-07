//! Positive orders checked against independent parameter integration.
use oneloop::{B0, OneLoopExpressions};
use symbolica::prelude::RealLike;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::{Complex, Float, SingleFloat},
};

#[test]
fn b0_epsilon_matches_parameter_integrals_and_mass_exchange() {
    std::thread::Builder::new().stack_size(128 * 1024 * 1024).spawn(|| {
        let args = ["eps_p", "eps_m0", "eps_m1", "eps_scale"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let call = B0().call((1, &args[0], &args[1], &args[2], &args[3]));
        let hook = Atom::evaluator_multiple(&[call.as_view()], &args)
            .direct_translation(true).build().unwrap();
        let map = OneLoopExpressions::new().evaluator(&[call], &args).unwrap();
        for (route, exact) in [("hook", hook), ("map", map)] {
            for bits in [53, 256] {
                let prototype = Float::new(bits);
                let mut high = exact.clone().map_coeff_with_prec(&|z| Complex::new(
                    prototype.from_rational(&z.re), prototype.from_rational(&z.im)), bits);
                let mut machine = exact.clone().map_coeff(&|z|Complex::new(z.re.to_f64(),z.im.to_f64()));
                for row in include_str!("b0_epsilon.tsv").lines().filter(|s|!s.starts_with('#')) {
                    let v: Vec<_> = row.split_whitespace().map(|s|Float::parse(s,Some(bits)).unwrap()).collect();
                    let z = || Float::new(bits);
                    let expected = Complex::new(v[6].clone(),v[7].clone());
                    for swap in [false,true] {
                        let (a,b) = if swap {(3,1)} else {(1,3)};
                        let inputs = [Complex::new(v[0].clone(),z()),Complex::new(v[a].clone(),v[a+1].clone()),
                            Complex::new(v[b].clone(),v[b+1].clone()),Complex::new(v[5].clone(),z())];
                        let mut result = [Complex::new(z(),z())];
                        high.evaluate(&inputs,&mut result);
                        let delta = &result[0] - &expected;
                        let err = delta.re.to_f64().hypot(delta.im.to_f64());
                        let tolerance = if bits==53 {2e-8} else {1e-65};
                        assert!(err<tolerance,"{route}, {bits} bits, swap={swap}, {row}: {} vs {expected}, error {err}",result[0]);
                        if bits==53 {
                            let mut out = [Complex::new(0.,0.)];
                            machine.evaluate(&inputs.map(|z|Complex::new(z.re.to_f64(),z.im.to_f64())),&mut out);
                            let err=(out[0].re-expected.re.to_f64()).hypot(out[0].im-expected.im.to_f64());
                            assert!(err<2e-8,"{route}, f64, swap={swap}, {row}: {:?}, error {err}",out[0]);
                        }
                    }
                }
            }
        }
        // Compiled transparent definitions use the same coefficient and branch convention.
        let call = B0().call((1, &args[0], &args[1], &args[2], &args[3]));
        let mut jit = OneLoopExpressions::new().jit_evaluator(&[call], &args).unwrap();
        for row in include_str!("b0_epsilon.tsv").lines().filter(|s|!s.starts_with('#')) {
            let v: Vec<f64> = row.split_whitespace().map(|s|s.parse().unwrap()).collect();
            let inputs = [Complex::new(v[0],0.),Complex::new(v[1],v[2]),Complex::new(v[3],v[4]),Complex::new(v[5],0.)];
            let mut output = [Complex::new(0.,0.)];
            jit.evaluate(&inputs, &mut output).unwrap();
            let err=(output[0].re-v[6]).hypot(output[0].im-v[7]);
            assert!(err<2e-8,"JIT, {row}: {:?}, error {err}",output[0]);
        }
        // Construction with exact degenerate arguments must also work.
        for (p,a,b) in [(0,0,0),(0,2,2),(4,1,1),(2,0,2)] {
            let expression = oneloop::b0_epsilon(&Atom::num(p), &Atom::num(a), &Atom::num(b), &Atom::num(1));
            let mut evaluator = OneLoopExpressions::new().evaluator(&[expression], &[]).unwrap().map_coeff(&|z| Complex::new(z.re.to_f64(), z.im.to_f64()));
            let mut output = [Complex::new(0., 0.)];
            evaluator.evaluate(&[], &mut output);
            assert!(output[0].re.is_finite() && output[0].im.is_finite());
        }
    }).unwrap().join().unwrap();
}
