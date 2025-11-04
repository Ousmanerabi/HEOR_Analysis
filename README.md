# Lung Transplantation – HEOR Portfolio Script (Synthetic Data)

This R script demonstrates a complete **Health Economics and Outcomes Research (HEOR)** workflow applied to a **synthetic cohort** of lung and heart–lung transplant patients.  
It was developed for educational and portfolio purposes to illustrate cost–utility analysis using parametric survival modeling.

---

## 🎯 Objective

The script reproduces the main analytical components of a real-world study performed at *Foch Hospital (France)*, while using **fully simulated data**.  
It aims to evaluate the **cost-effectiveness** of lung transplantation compared with medical management.

---

## ⚙️ Workflow Overview

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

## 🧠 Key Methods

- **Survival model:** Weibull (scale ≈ 1670 days, shape ≈ 1.36)  
- **Perspective:** French National Health Insurance  
- **Outcome:** Cost per QALY gained  
- **Sensitivity:** One-way and PSA (1,000 simulations)

---

## 🧩 How to Run

```r
# 1. Install dependencies (first time)
install.packages(c("tidyverse", "survival", "flexsurv", "MASS"))

# 2. Source the script
source("lung_transplant_HEOR_portfolio.R")
