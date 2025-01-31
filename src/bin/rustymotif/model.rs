/// A Beta-Bernoulli model
#[derive(Clone, Debug)]
pub struct BetaBernoulliModel {
    pub alpha: f64,
    pub beta: f64,
}

impl BetaBernoulliModel {
    /// Create a new Beta-Bernoulli model with alpha=1.0, beta=1.0
    pub fn new() -> Self {
        Self {
            alpha: 1.0,
            beta: 1.0,
        }
    }

    pub fn new_with_params(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    /// Update parameters with `n_meth` successes and `n_nonmeth` failures
    pub fn update(&mut self, n_meth: usize, n_nonmeth: usize) {
        self.alpha += n_meth as f64;
        self.beta += n_nonmeth as f64;
    }

    /// Compute mean of posterior
    pub fn mean(&self) -> f64 {
        self.alpha / (self.alpha + self.beta)
    }

    /// Compute standard deviation of posterior
    pub fn standard_deviation(&self) -> f64 {
        let numerator = self.alpha * self.beta;
        let denominator = (self.alpha + self.beta).powf(2.0) * (self.alpha + self.beta + 1.0);
        (numerator / denominator).sqrt()
    }
}
