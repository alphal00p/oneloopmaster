use symbolica::{domains::float::Real, evaluate::ExpressionEvaluator};

// Upstream Symbolica provides a scalar interpreter. These row oracles preserve
// mathematical and branch checks, but do not exercise the fork's removed batch
// interpreter or make claims about its instruction order or workspace reuse.
pub trait EvaluateRows<T> {
    fn evaluate_rows(&mut self, input: &[T], output: &mut [T], rows: usize) -> Result<(), String>;
}

impl<T: Real> EvaluateRows<T> for ExpressionEvaluator<T> {
    fn evaluate_rows(&mut self, input: &[T], output: &mut [T], rows: usize) -> Result<(), String> {
        if rows == 0 {
            return if input.is_empty() && output.is_empty() {
                Ok(())
            } else {
                Err("zero rows require empty input and output buffers".into())
            };
        }
        if input.len() % rows != 0 || output.len() % rows != 0 {
            return Err("incomplete evaluator row".into());
        }
        let inputs = input.len() / rows;
        let outputs = output.len() / rows;
        for row in 0..rows {
            self.try_evaluate(
                &input[row * inputs..(row + 1) * inputs],
                &mut output[row * outputs..(row + 1) * outputs],
            )
            .map_err(|error| error.to_string())?;
        }
        Ok(())
    }
}
