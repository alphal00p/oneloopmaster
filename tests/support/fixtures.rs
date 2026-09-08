//! Shared parser for the original-Fortran scalar fixture format.
use oneloop::ScalarIntegral;
use symbolica::domains::float::Complex;

pub struct Fixture {
    pub family: ScalarIntegral,
    pub name: String,
    pub args: Vec<Complex<f64>>,
    pub expected: [Complex<f64>; 3],
    pub normalization: f64,
}

pub fn parse(source: &str) -> Vec<Fixture> {
    let mut name = String::new();
    let mut rows = Vec::new();
    for (line_number, line) in source.lines().enumerate() {
        if let Some(rest) = line.strip_prefix("# name ") {
            name = rest.split_whitespace().next().unwrap().into();
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let values: Vec<f64> = line
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        let (family, momenta, masses, power) = match values[0] {
            1. => (ScalarIntegral::A0, 0, 1, -1),
            2. => (ScalarIntegral::B0, 1, 2, 0),
            -2. => (ScalarIntegral::DB0, 1, 2, 1),
            3. => (ScalarIntegral::C0, 3, 3, 1),
            4. => (ScalarIntegral::D0, 6, 4, 2),
            _ => panic!("invalid scalar fixture family"),
        };
        let mass_start = 2 + momenta;
        let reference_start = mass_start + 2 * masses;
        assert_eq!(values.len(), reference_start + 6);
        let mut args: Vec<_> = values[2..mass_start]
            .iter()
            .map(|&v| Complex::new(v, 0.))
            .collect();
        args.extend(
            values[mass_start..reference_start]
                .chunks_exact(2)
                .map(|v| Complex::new(v[0], v[1])),
        );
        args.push(Complex::new(values[1], 0.));
        // Match the Fortran survey: independent scale variation must not change
        // the absolute-error allowance for identical physical kinematics.
        let scale = values[2..reference_start]
            .iter()
            .fold(0.0_f64, |scale, value| scale.max(value.abs()));
        let scale = if scale == 0. { 1. } else { scale };
        rows.push(Fixture {
            family,
            name: if name.is_empty() {
                format!("line_{}", line_number + 1)
            } else {
                std::mem::take(&mut name)
            },
            args,
            expected: core::array::from_fn(|i| {
                Complex::new(
                    values[reference_start + 2 * i],
                    values[reference_start + 2 * i + 1],
                )
            }),
            normalization: scale.powi(power),
        });
    }
    rows
}

impl Fixture {
    #[allow(dead_code)] // Some clients aggregate failures instead of stopping at the first row.
    pub fn check(&self, output: &[Complex<f64>]) {
        let failures = self.failures(output);
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }

    pub fn failures(&self, output: &[Complex<f64>]) -> Vec<String> {
        assert_eq!(output.len(), 3);
        let mut failures = Vec::new();
        for (coefficient, (actual, expected)) in output.iter().zip(self.expected).enumerate() {
            let error =
                (actual.re - expected.re).hypot(actual.im - expected.im) * self.normalization;
            let tolerance = 2e-10 + 2e-8 * expected.re.hypot(expected.im) * self.normalization;
            if !actual.re.is_finite()
                || !actual.im.is_finite()
                || !error.is_finite()
                || error > tolerance
            {
                failures.push(format!("{} {} coefficient {coefficient}: {actual:?}, expected {expected:?}, normalized error {error}", self.family.name(), self.name));
            }
        }
        failures
    }
}
