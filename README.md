# Lung Transplantation – Cost-Utility-Analysis

This project demonstrates a full end-to-end Health Economics and Outcomes Research (HEOR) workflow using parametric survival modeling and cost–utility analysis in R.

It replicates the analytical framework of a real-world lung transplantation study conducted at Foch Hospital (France), using a fully simulated cohort for reproducibility.

---

##  Objective

To evaluate the cost-effectiveness of lung transplantation compared to medical management from the payer perspective (French National Health Insurance), using survival modeling and QALY estimation.

---

##  Workflow Overview

1. **Data Simulation**  
   - Generates a synthetic cohort of 631 patients.  
   - Includes demographics, clinical indicators, transplant status, and quality-of-life (QoL) data.

2. **Survival Analysis**  
   - Fits a **parametric Weibull model** to post-transplant survival data.  
   - Estimates median and mean survival time and 1–5 year survival probabilities.

3. **Quality of Life (QoL)**  
   - Simulates NHP-based dimensions and maps them to **EQ-5D utilities** (pre- and post-transplant).  
   - Calculates incremental **QALYs**.

4. **Cost–Utility Analysis**  
   - Computes the **Incremental Cost-Effectiveness Ratio (ICER)** from the payer perspective (Assurance Maladie).  
   - Includes deterministic and probabilistic sensitivity analyses.

---

## Key Methods

- **Survival model:** Weibull (scale ≈ 1670 days, shape ≈ 1.36)  
- **Perspective:** French National Health Insurance  
- **Outcome:** Cost per QALY gained  
- **Sensitivity:** One-way and PSA (1,000 simulations)

---

## Key Outputs

- Mean and median survival estimates

- Incremental QALYs

- ICER (Cost per QALY gained)

- Tornado diagram

- PSA distributions


## What This Project Demonstrates

- Advanced parametric survival modeling in R (flexsurv)

- Health technology assessment (HTA) methodology

- ICER computation and payer-perspective modeling

- Deterministic and probabilistic uncertainty quantification

- Reproducible analytical pipelines

- Clear separation between data generation, modeling, and economic evaluation steps

## Technical Stack

R packages used:

- tidyverse

- survival

- flexsurv

- MASS


##  How to Run

```r
# 1. Install dependencies (first time)
install.packages(c("tidyverse", "survival", "flexsurv", "MASS"))

# 2. Source the script
source("lung_transplant_HEOR_portfolio.R")
