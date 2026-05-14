
This README file was generated on **2026-05-14** by **Brock D. Sherlock**

---

# GENERAL INFORMATION

## Title of Dataset:
Code and synthetic data for *Does Timing Matter? Exploring the Effects of Measurement Error on Models*

---

## Author Information  
Name: Brock D. Sherlock  
ORCID: 0000-0001-9072-5176 
Institution: University of New South Wales (UNSW)  
Address: School of Mathematics & Statistics, UNSW, Sydney NSW 2052, Australia  
Email: brock.sherlock@unsw.edu.au  

---

## Author/Alternate Contact Information  
Name: Adelle C. F. Coster  
ORCID: 0000-0002-5572-6832 
Institution: University of New South Wales (UNSW)  
Address: School of Mathematics & Statistics, UNSW, Sydney NSW 2052, Australia  
Email: A.Coster@unsw.edu.au  

---


* Date of data collection: 
2025 (year synthetic data and analyses produced)
* Geographic location of data collection: 
Sydney, NSW, Australia
* Information about funding sources that supported the collection of the data: 
Australian Research Council (ARC) Discovery Projects funding scheme (DP210100255)


# SHARING/ACCESS INFORMATION

* Links to publications that cite or use the data
Sherlock et al. (2026), *Does Timing Matter? Exploring the Effects of Measurement Error on Models*, Bulletin of Mathematical Biology

* Links to other publicly accessible locations of the data: 
GitHub repository:  
https://github.com/brocksherlockmaths/Does-Timing-Matter


# DATA & FILE OVERVIEW

## File List

### Core Scripts

- `linearRegression_confidenceInterval_plot.m`  
  Demonstrates effects of classical and Berkson measurement error in linear regression. Produces fitted models and confidence intervals under controlled and non-controlled experimental designs.

- `caseStudy_oscillation_ci.m`  
  Generates oscillatory datasets and compares fitted models with confidence intervals under multiple measurement error scenarios.

- `caseStudy_oscillation_density.m`  
  Performs large-scale simulations (10,000 repeats) to estimate distributions of parameter estimates for oscillatory systems.

- `parasitemia_caseStudy_ci.m`  
  Simulates parasite dynamics using a nonlinear biological model and investigates the impact of observation and time measurement error on parameter estimation.

- `caseStudy_tumour_ci.m`  
  Fits Gompertz tumour growth models to real and synthetic data and compares confidence intervals under different timing error scenarios.

- `caseStudy_tumour.m`  
  Performs full tumour growth simulations including:
  - synthetic data generation  
  - correlated time error modelling  
  - repeated simulations (10,000 replicates)  
  - parameter inference and uncertainty analysis  

---

### Helper Functions

- `generate_correlatedTimes.m`  
  Generates correlated timing errors using cumulative random delays.  
  Supports:
  - lognormal timing errors (default)
  - normal timing errors  
  Produces:
  - ordered measurement times (fixed sampling sequence)
  - randomised measurement times (random sampling order)

- `normalGenerator_noTimeError.m`  
  Generates synthetic data with observational noise only (no timing error).  
  Replicates model output and adds Gaussian noise independently at each time point.

- `normalGenerator_normalTimeError.m`  
  Generates synthetic data with additive Gaussian timing error.  
  Time points are perturbed before evaluating the model, and observational noise is added afterward.

---

# METHODOLOGICAL INFORMATION

## Description of methods used for collection/generation of data

All data are synthetically generated using mathematical models with controlled measurement error.

Models include:
- Linear regression models  
- Oscillatory cosine models  
- Parasite growth models  
- Tumour growth models (Gompertz equation)  
- GLUT4 translocation compartmental model  

Measurement error is introduced via:
- Classical error: \(W = X + U\)  
- Berkson error: \(X = W + U\)  

Simulation studies include repeated experiments (up to 10,000 datasets) to evaluate parameter bias and variability.

---

## Methods for processing the data

- Data generated from model equations  
- Measurement error applied using specified distributions  
- Parameter inference performed using:
  - Least squares
  - Maximum likelihood estimation
  - Bayesian inference (Sequential Monte Carlo)  
- Outputs include parameter estimates, confidence intervals, and plots  

---

## Instrument- or software-specific information needed

- MATLAB (standard environment)  
- Core MATLAB functions used:
  - `fit`
  - `confint`
  - `ksdensity`  

Bayesian analyses implemented using Sequential Monte Carlo methods in MATLAB scripts.

## Specialized formats or other abbreviations used

- MLE: Maximum Likelihood Estimation  
- SMC: Sequential Monte Carlo  
- MAP: Maximum A Posteriori  
- CI: Confidence Interval  

---

