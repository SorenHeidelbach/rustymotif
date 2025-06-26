use statrs::function::gamma::ln_gamma;
use anyhow::Result;
use crate::pileup::PileupRecord;

use log::{debug, warn};
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


#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Beta {
    pub alpha: f64,
    pub beta: f64,
}

impl Beta {
    pub fn new(alpha: f64, beta: f64) -> Self {
        Self { alpha, beta }
    }

    pub fn from_data(records: &Vec<&PileupRecord>) -> Self {
        let mut mean = 0.0;
        let mut sum_sq = 0.0;
        let mut count = 0.0;
        for record in records {
            if record.n_valid_cov == 0 {
                continue;
            }
            let degree = record.n_mod as f64 / record.n_valid_cov as f64;
            mean += degree;
            sum_sq += degree * degree;
            count += 1.0;
        }
        mean /= count;
        let variance = sum_sq / count - mean * mean;
        let common = mean * (1.0 - mean) / variance - 1.0;
        let alpha = mean * common;
        let beta = (1.0 - mean) * common;
        Beta::new(alpha, beta)
    }

    pub fn mean(&self) -> f64 {
        self.alpha / (self.alpha + self.beta)
    }

    pub fn variance(&self) -> f64 {
        let numerator = self.alpha * self.beta;
        let denominator = (self.alpha + self.beta).powf(2.0) * (self.alpha + self.beta + 1.0);
        numerator / denominator
    }

    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }

    pub fn log_beta(&self) -> f64 {
        ln_gamma(self.alpha) + ln_gamma(self.beta) - ln_gamma(self.alpha + self.beta)
    }

    pub fn log_pdf(&self, x: f64) -> f64 {
        if x < 0.0 || x > 1.0 {
            return f64::NEG_INFINITY;
        }
        (self.alpha - 1.0) * x.ln() + (self.beta - 1.0) * (1.0 - x).ln() - self.log_beta()
    }
    pub fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 || x > 1.0 {
            return 0.0;
        }
        self.log_pdf(x).exp()
    }
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BetaMixture {
    pub true_model: Beta,
    pub false_model: Beta,
    pub pi: f64,
}
impl BetaMixture {
    pub fn new(true_model: Beta, false_model: Beta, pi: f64) -> Self {
        Self {
            true_model,
            false_model,
            pi,
        }
    }

    pub fn log_pdf(&self, x: f64) -> f64 {
        let log_p_true = self.true_model.log_pdf(x);
        let log_p_false = self.false_model.log_pdf(x);
        let log_term_true = self.pi.ln() + log_p_true;
        let log_term_false = (1.0 - self.pi).ln() + log_p_false;
        log_sum_exp(log_term_true, log_term_false)
    }
}
impl BetaMixture {
    pub fn mean(&self) -> f64 {
        self.true_model.mean() * self.pi + self.false_model.mean() * (1.0 - self.pi)
    }

    pub fn variance(&self) -> f64 {
        let true_var = self.true_model.variance();
        let false_var = self.false_model.variance();
        let true_mean = self.true_model.mean();
        let false_mean = self.false_model.mean();
        let pi = self.pi;
        pi * (true_var + (true_mean - self.mean()).powi(2)) +
            (1.0 - pi) * (false_var + (false_mean - self.mean()).powi(2))
    }

    pub fn log_likelihoods(
        &self,
        records: &Vec<&PileupRecord>,
    ) -> (f64, f64) {
        let mut log_lik_true = 0.0;
        let mut log_lik_false = 0.0;
        let epsilon = 1e-6;
        for record in records {
            if record.n_valid_cov == 0 {
                continue;
            }
            let raw_degree = (record.n_mod as f64) / (record.n_valid_cov as f64);
            let degree = raw_degree.max(epsilon).min(1.0 - epsilon);
            log_lik_true += self.true_model.log_pdf(degree);
            log_lik_false += self.false_model.log_pdf(degree);
        }
        (log_lik_true, log_lik_false)
    }

    // Given values and corresponding weights, this function computes the weighted
    // mean and variance and then returns a Beta distribution estimated using the method-of-moments.
    pub fn beta_weighted(
        &self,
        values: &Vec<f64>,
        weights: &Vec<f64>,
        prior_alpha: f64,
        prior_beta: f64,
    ) -> Option<Beta> {
        let mut weight_sum = 0.0;
        let mut weighted_sum = 0.0;
        for (&value, &weight) in values.iter().zip(weights.iter()) {
            weight_sum += weight;
            weighted_sum += weight * value;
        }
        if weight_sum == 0.0 {
            debug!("No weights provided, cannot fit Beta distribution.");
            return None;
        }
        let mean = weighted_sum / weight_sum;

        let mut weighted_var = 0.0;
        for (&value, &weight) in values.iter().zip(weights.iter()) {
            weighted_var += weight * (value - mean).powi(2);
        }
        weighted_var /= weight_sum;
        if weighted_var < 1e-12 {
            debug!("Variance is too small, cannot fit Beta distribution.");
            return None;
        }
        
        let common = mean * (1.0 - mean) / weighted_var - 1.0;
        if common <= 0.0 {
            debug!("Common term is non-positive, cannot fit Beta distribution.");
            return None;
        }

        // Method-of-moments estimates
        let alpha_est = mean * common;
        let beta_est = (1.0 - mean) * common;

        // Incorporate the prior by adding pseudo-counts.
        let new_alpha = prior_alpha + alpha_est;
        let new_beta = prior_beta + beta_est;
        Some(Beta::new(new_alpha, new_beta))
    }

    pub fn fit_em_fixed_false(
        &mut self,
        records: &Vec<&PileupRecord>,
        false_penalty: f64,
        max_iter: usize,
        tol: f64,
    ) -> Result<f64> {
        let mut log_lik_true_old:f64;
        let mut log_lik_true:f64 = 0.0;
        let epsilon = 1e-12;
    
        for _ in 0..max_iter {
            let mut responsibilities = Vec::new();
            log_lik_true = 0.0;
    
            // E-step
            for record in records {
                if record.n_valid_cov == 0 {
                    continue;
                }
                let raw_degree = (record.n_mod as f64) / (record.n_valid_cov as f64);
                let degree = raw_degree.max(epsilon).min(1.0 - epsilon);

                let log_p_true = match self.true_model.log_pdf(degree) {
                    log_p if log_p.is_finite() => log_p,
                    _ => {
                        warn!("Invalid log_pdf for true model at degree: {}", degree);
                        continue;
                    }
                };


                let log_p_false = match self.false_model.log_pdf(degree) {
                    log_p if log_p.is_finite() => log_p,
                    _ => {
                        warn!("Invalid log_pdf for false model at degree: {}", degree);
                        continue;
                    }
                };
                
                let log_term_true = self.pi.ln() + log_p_true;
                let penalty = 1.0 / (1.0 + ((degree - 0.25) * 10.0).exp());
                let log_penalty = penalty.ln();
                let log_term_false = (1.0 - self.pi).ln() + log_p_false + log_penalty;
                let log_mix = log_sum_exp(log_term_true, log_term_false);
                let log_responsibility = log_term_true - log_mix;
                let responsibility = log_responsibility.exp();
                responsibilities.push((degree, responsibility));

                log_lik_true += log_mix;
            }
    
            // M-step:
            if responsibilities.is_empty() {
                warn!("No valid records found for fitting the model.");
                return Err(anyhow::anyhow!("No valid records found for fitting the model."));
            }
            let new_pi = responsibilities.iter().map(|&(_, r)| r).sum::<f64>()
                / (responsibilities.len() as f64);
    
            let (mut true_degrees, mut true_weights) = (Vec::new(), Vec::new());
            for &(degree, r) in &responsibilities {
                true_degrees.push(degree);
                true_weights.push(r);
            }
            self.true_model = match self.beta_weighted(&true_degrees, &true_weights, 0.0, 0.0) {
                Some(model) => model,
                None => {
                    warn!("Failed to fit true model with weighted degrees");
                    debug!("True degrees: {:?}", true_degrees);
                    debug!("True weights: {:?}", true_weights);
                    debug!("False model: α = {:.6}, β = {:.6}", self.false_model.alpha, self.false_model.beta);
                    debug!("π: {:.6}", self.pi);
                    debug!("Log likelihood: {:.6}", log_lik_true);
                    debug!("Responsibilities: {:?}", responsibilities);
                    debug!("Responsibilities sum: {:?}", responsibilities.iter().map(|&(_, r)| r).sum::<f64>());
                    debug!("Responsibilities length: {:?}", responsibilities.len());
                    
                    return Err(anyhow::anyhow!("Failed to fit true model"));
                }
            };
            self.pi = new_pi;
            log_lik_true_old = log_lik_true;
            
            // Termination condition
            if (log_lik_true - log_lik_true_old).abs() < tol {
                // debug!("Final log_likelihood: {:.3}, π: {:.3}, α = {:.3}, β = {:.3}, mean: {:.2}", log_lik_true, self.pi, self.true_model.alpha, self.true_model.beta, self.true_model.mean());
                break;
            }
        }
    
        Ok(log_lik_true)
    }

}




fn weighted_beta_fit_degrees_with_prior(degrees: &[f64], weights: &[f64], prior_alpha: f64, prior_beta: f64) -> Option<Beta> {
    let mut weight_sum = 0.0;
    let mut weighted_sum = 0.0;
    for (&d, &w) in degrees.iter().zip(weights.iter()) {
        weight_sum += w;
        weighted_sum += w * d;
    }
    if weight_sum == 0.0 {
        return None;
    }
    let mean = weighted_sum / weight_sum;
    let mut weighted_var = 0.0;
    for (&d, &w) in degrees.iter().zip(weights.iter()) {
        weighted_var += w * (d - mean).powi(2);
    }
    weighted_var /= weight_sum;
    // Avoid division by zero or very low variance.
    if weighted_var < 1e-12 {
        return None;
    }
    let common = mean * (1.0 - mean) / weighted_var - 1.0;
    if common <= 0.0 {
        return None;
    }
    // Method-of-moments estimates (based solely on the observed data)
    let alpha_est = mean * common;
    let beta_est = (1.0 - mean) * common;
    // Incorporate the prior by adding pseudo-counts.
    let new_alpha = prior_alpha + alpha_est;
    let new_beta = prior_beta + beta_est;
    Some(Beta::new(new_alpha, new_beta))
}

fn log_sum_exp(a: f64, b: f64) -> f64 {
    match (a, b) {
        (x, y) if x.is_infinite() && x.is_sign_negative() => y,
        (x, y) if y.is_infinite() && y.is_sign_negative() => x,
        (x, y) => {
            let max_val = x.max(y);
            max_val + ((x - max_val).exp() + (y - max_val).exp()).ln()
        }
    }
}
// --- EM Algorithm for Fitting Only the True Beta Component ---
//
// In this version, the false component is fixed (known) and we only update the true component and the mixing weight π.
// We compute the methylation degree for each record (clamped to avoid boundary issues).
// In the E-step, we compute the log likelihoods for the true and false models using the log_pdf method,
// then combine them using the log-sum-exp trick. We then add a penalty term to the false component so that
// observations that are well explained by the false model contribute even less to the responsibility for the true model.
// This penalty term is added in log-space.
pub fn fit_true_beta_model_em(
    records: &Vec<&PileupRecord>,
    false_model: &Beta,
) -> (Beta, f64, f64, f64) {
    debug!("Fitting true beta model using EM algorithm");
    debug!("False model: α = {:.6}, β = {:.6}", false_model.alpha, false_model.beta);
    

    let mut true_model = Beta::new(5.0, 1.0);
    let mut pi = 0.05_f64; 

    let max_iter = 100;
    let tol = 1e-6;
    let mut log_lik_old = f64::NEG_INFINITY;
    let mut log_lik = 0.0;
    let mut log_lik_false = 0.0;

    // Precompute the total log likelihood under the false model for reference.
    for record in records {
        if record.n_valid_cov == 0 {
            continue;
        }
        let raw_degree = (record.n_mod as f64) / (record.n_valid_cov as f64);
        let degree = raw_degree.max(1e-6).min(1.0 - 1e-6);
        log_lik_false += false_model.log_pdf(degree);
    }

    let epsilon = 1e-6; // For clamping degree values.
    let false_penalty = 2.0_f64; // Factor >1 to penalize observations explained by the false model.

    for iter in 0..max_iter {
        let mut responsibilities = Vec::new();
        log_lik = 0.0;

        // E-step: compute responsibilities in log-space.
        for record in records {
            if record.n_valid_cov == 0 {
                continue;
            }
            let raw_degree = (record.n_mod as f64) / (record.n_valid_cov as f64);
            let degree = raw_degree.max(epsilon).min(1.0 - epsilon);
            let log_p_true = true_model.log_pdf(degree);
            let log_p_false = false_model.log_pdf(degree);


            // Add penalty to the false component to downweight it.
            let log_penalty = false_penalty.ln(); // positive constant
            let log_term_true = pi.ln() + log_p_true;
            let log_term_false = (1.0 - pi).ln() + log_p_false - log_penalty;
            let log_mix = log_sum_exp(log_term_true, log_term_false);
            let mix = log_mix.exp(); // Stable mixture probability.
            let log_r = log_term_true - log_mix;
            let r = log_r.exp();
            responsibilities.push((degree, r));
            log_lik += log_mix;
        }

        // M-step:
        let new_pi = responsibilities.iter().map(|&(_, r)| r).sum::<f64>()
            / (responsibilities.len() as f64);

        let (mut true_degrees, mut true_weights) = (Vec::new(), Vec::new());
        for &(degree, r) in &responsibilities {
            true_degrees.push(degree);
            true_weights.push(r);
        }

        // Update the true model using a weighted Beta fit that incorporates a prior (here, Beta(5,1)).
        if let Some(new_true_model) = weighted_beta_fit_degrees_with_prior(&true_degrees, &true_weights, 5.0, 1.0) {
            true_model = new_true_model;
        }

        //debug!("Iteration: {}, log_likelihood: {:.6}, pi: {:.6}, true_model: α = {:.6}, β = {:.6}",
        //       iter, log_lik, new_pi, true_model.alpha, true_model.beta);

        if (log_lik - log_lik_old).abs() < tol {
            debug!("Final log_likelihood: {:.3}, π: {:.3}, α = {:.3}, β = {:.3}, mean: {:.2}", log_lik, new_pi, true_model.alpha, true_model.beta, true_model.mean());
            pi = new_pi;
            log_lik_old = log_lik;
            break;
        }
        pi = new_pi;
        log_lik_old = log_lik;
        //debug!("Iteration: {}, log_likelihood: {:.6}, pi: {:.6}, α = {:.6}, β = {:.6}",
        //       iter, log_lik, pi, true_model.alpha, true_model.beta);
    }

    (true_model, pi, log_lik, log_lik_false)
}