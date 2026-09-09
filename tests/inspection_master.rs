//! Expression-first inspection, including composite inputs and native attributes.
use oneloop::{ExpressionOptions, ScalarIntegral, get_expression, get_expression_with_options};
use symbolica::atom::{Atom, AtomCore};

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn master_calls_accept_composite_kinematics_and_native_assumptions() {
    on_stack(|| {
        let s = symbolica::symbol!("master_inspection_s"; Real).to_atom();
        let m = symbolica::symbol!("master_inspection_m"; Positive).to_atom();
        let u = symbolica::symbol!("master_inspection_u"; Positive).to_atom();
        let input = [s.pow(2) - 3, &m + 4, m.pow(2), &u + 1];
        let call = oneloop::B0().call(input.as_slice());
        let inspected = get_expression(&call).unwrap();
        assert_eq!(
            inspected,
            ScalarIntegral::B0.get_expression(&input).unwrap()
        );
        assert_eq!(get_expression(call.as_view()).unwrap(), inspected);
        assert_eq!(
            oneloop::master_arguments(&call).unwrap(),
            (ScalarIntegral::B0, input.to_vec())
        );
        let alias = symbolica::symbol!("oneloop::B0").call(input.as_slice());
        assert_eq!(get_expression(alias).unwrap(), inspected);
    });
}

#[test]
fn master_calls_cover_all_families_and_validate_shape() {
    on_stack(|| {
        for (family, master, momenta) in [
            (ScalarIntegral::A0, oneloop::A0(), 0),
            (ScalarIntegral::B0, oneloop::B0(), 1),
            (ScalarIntegral::DB0, oneloop::dB0(), 1),
            (ScalarIntegral::C0, oneloop::C0(), 3),
            (ScalarIntegral::D0, oneloop::D0(), 6),
        ] {
            let input = (0..family.arity())
                .map(|i| Atom::num(i32::from(i >= momenta)))
                .collect::<Vec<_>>();
            assert_eq!(
                get_expression(master.call(&input)).unwrap(),
                family.get_expression(&input).unwrap()
            );
        }
        let x = symbolica::symbol!("master_inspection_invalid").to_atom();
        assert!(get_expression(&x).unwrap_err().contains("master call"));
        assert!(get_expression(&x + 1).is_err());
        assert!(
            get_expression(symbolica::symbol!("unrelated::B0").call((0, 1, 1, 1)))
                .unwrap_err()
                .contains("not a supported")
        );
        assert!(
            get_expression(oneloop::B0().call((0, 1, 1)))
                .unwrap_err()
                .contains("expects 4")
        );
        assert!(
            get_expression(oneloop::B0().call((0, 0, 1, 1, 1)))
                .unwrap_err()
                .contains("no Laurent tag")
        );
        assert!(
            get_expression_with_options(
                oneloop::A0().call((&x, 1)),
                ExpressionOptions {
                    max_nodes: 0,
                    max_depth: 512
                }
            )
            .is_err()
        );
    });
}
