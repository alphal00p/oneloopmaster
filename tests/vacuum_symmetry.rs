//! Independent analytic vacuum checks; unstable Fortran orders are not targets.
use oneloop::OneLoopExpressions;
use symbolica::prelude::*;

#[test]
fn vacuum_triangle_preserves_mass_permutations_and_lower_lips() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let x = ["vac_m1", "vac_m2", "vac_m3", "vac_mu"].map(|s| symbol!(s).to_atom());
            let zero = Atom::num(0);
            let context = OneLoopExpressions::new();
            let series = context.c0([&zero; 3], [&x[0], &x[1], &x[2]], &x[3]);
            let exact = context.evaluator(series.coefficients(), &x).unwrap();
            let permutations = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]];
            for bits in [128, 256] {
                let mut evaluator = exact.clone().map_coeff_with_prec(&|c| Complex::new(c.re.to_multi_prec_float(bits), c.im.to_multi_prec_float(bits)), bits);
                for real in [[-1.,2.,3.],[1.,2.,3.],[-3.,-2.,-1.]] {
                    for width in [0., 0.1, 1e-12] {
                        let masses = real.map(|r| Complex::new(Float::with_val(bits,r),Float::with_val(bits,-width)));
                        // Partial fractions of the three propagators give
                        // -sum m_i log(m_i-i0)/prod_{j!=i}(m_i-m_j).
                        let mut expected = Complex::new(Float::new(bits),Float::new(bits));
                        for i in 0..3 {
                            let mut logarithm = masses[i].clone().log();
                            if width == 0. && real[i] < 0. {
                                logarithm = Complex::new(Float::with_val(bits,-real[i]).log(),-Float::with_val(bits,symbolica::domains::backend::float::Constant::Pi));
                            }
                            let denominator = (&masses[i]-&masses[(i+1)%3])*(&masses[i]-&masses[(i+2)%3]);
                            expected -= masses[i].clone()*logarithm/denominator;
                        }
                        for permutation in permutations {
                            for mu_squared in [0.01,1.,10000.] {
                                let input = [masses[permutation[0]].clone(),masses[permutation[1]].clone(),masses[permutation[2]].clone(),Complex::new(Float::with_val(bits,mu_squared),Float::new(bits))];
                                let mut output = core::array::from_fn::<_,3,_>(|_|Complex::new(Float::new(bits),Float::new(bits)));
                                evaluator.evaluate(&input,&mut output);
                                let delta = output[0].clone()-&expected;
                                let error = delta.re.to_f64().hypot(delta.im.to_f64());
                                let tolerance = if bits==128 {1e-28} else {1e-60};
                                assert!(error.is_finite() && error<tolerance,"masses={real:?},width={width},permutation={permutation:?},mu²={mu_squared},bits={bits},actual={:?},expected={expected:?},error={error}",output[0]);
                                assert!(output[1].is_zero() && output[2].is_zero());
                            }
                        }
                    }
                }
            }
        })
        .unwrap().join().unwrap();
}
