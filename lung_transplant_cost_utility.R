############################################################
# Lung / Heart–Lung Transplantation – Cost-Utility-Analysis
# Author: Ousmane Diallo
# Description:
#   - Simulate a synthetic transplant cohort (no real data)
#   - Weibull survival modeling
#   - QoL utilities (EQ-5D) and QALY calculation
#   - Cost-utility analysis (ICER)
#   - One-way and probabilistic sensitivity analyses
############################################################

## 0. Packages ----
req <- c("tidyverse", "survival", "flexsurv", "MASS")
new <- setdiff(req, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(req, library, character.only = TRUE))

set.seed(2025)

############################################################
## 1. Simulate synthetic cohort (no real patient data) ----
############################################################

n_patients <- 631

# IDs
ID <- 1:n_patients

# Sex: 0 = female, 1 = male (≈50/50)
SEX <- rbinom(n_patients, size = 1, prob = 0.504)

# Age at listing: mean 38.8 (range approx 11–72)
AGE <- rnorm(n_patients, mean = 38.8, sd = 12)
AGE <- pmin(pmax(AGE, 11.6), 72.2)

# Diagnosis distribution (approximate percentages)
DIAGNOSIS <- sample(
  c("Cystic fibrosis", "COPD/emphysema", "Interstitial lung disease",
    "Other", "Pulmonary arterial hypertension"),
  size = n_patients,
  replace = TRUE,
  prob = c(0.506, 0.212, 0.152, 0.122, 0.008)
)

# Transplantation status: ~90.5% transplanted
TRANSPLANTED <- rbinom(n_patients, size = 1, prob = 0.905)

# Second transplant (~3% of transplanted)
SECOND_TX <- ifelse(TRANSPLANTED == 1,
                    rbinom(n_patients, size = 1, prob = 0.03),
                    0)

# Super-urgent status among transplanted (~10%)
SUPER_URGENT <- ifelse(TRANSPLANTED == 1,
                       rbinom(n_patients, size = 1, prob = 0.10),
                       0)

############################################################
## 2. Survival after transplantation (Weibull) ----
############################################################
# Only transplanted patients contribute to post-TP survival model

shape_tp <- 1.36     # shape parameter (k)
scale_tp <- 1670     # scale parameter (lambda) in days

n_tp <- sum(TRANSPLANTED == 1)

# True event times follow Weibull(k, lambda)
true_surv_tp <- rweibull(n_tp, shape = shape_tp, scale = scale_tp)

# Censoring times (e.g. administrative censoring at ~5 years)
cens_tp <- rexp(n_tp, rate = 1 / (5 * 365))

# Observed time and event indicator
time_obs_tp <- pmin(true_surv_tp, cens_tp)
event_tp    <- as.integer(true_surv_tp <= cens_tp)

# Initialize full vectors
DELAI <- rep(NA_real_, n_patients)   # time in days
EVENT <- rep(NA_integer_, n_patients)

# Fill transplanted patients
DELAI[TRANSPLANTED == 1] <- time_obs_tp
EVENT[TRANSPLANTED == 1] <- event_tp

# For non-transplanted, simulate time on waiting list (not used in main Weibull model)
DELAI[TRANSPLANTED == 0] <- rexp(sum(TRANSPLANTED == 0), rate = 1 / (2 * 365))
EVENT[TRANSPLANTED == 0] <- 0L   # censored

############################################################
## 3. QoL: NHP dimensions and EQ-5D utilities ----
############################################################

# We simulate QoL for a subset of transplanted patients (approx 302 like in the original study)
n_qol <- 302
tp_indices <- which(TRANSPLANTED == 1)
n_qol <- min(n_qol, length(tp_indices))
qol_ids <- sample(tp_indices, size = n_qol, replace = FALSE)

sim_nhp_dim <- function(n, mean_pct, sd_pct) {
  x <- rnorm(n, mean = mean_pct, sd = sd_pct)
  pmin(pmax(x, 0), 100)
}

# Initialize NHP PRE/POST
NHP_PAIN_PRE    <- NHP_PAIN_POST    <- rep(NA_real_, n_patients)
NHP_EMO_PRE     <- NHP_EMO_POST     <- rep(NA_real_, n_patients)
NHP_MOB_PRE     <- NHP_MOB_POST     <- rep(NA_real_, n_patients)
NHP_SLEEP_PRE   <- NHP_SLEEP_POST   <- rep(NA_real_, n_patients)
NHP_ENERGY_PRE  <- NHP_ENERGY_POST  <- rep(NA_real_, n_patients)
NHP_SOCIAL_PRE  <- NHP_SOCIAL_POST  <- rep(NA_real_, n_patients)

# Approximate means/SDs based on your tables (PRE)
NHP_PAIN_PRE[qol_ids]   <- sim_nhp_dim(n_qol, 15.5, 15)
NHP_EMO_PRE[qol_ids]    <- sim_nhp_dim(n_qol, 17.2, 18)
NHP_MOB_PRE[qol_ids]    <- sim_nhp_dim(n_qol, 35.3, 25)
NHP_SLEEP_PRE[qol_ids]  <- sim_nhp_dim(n_qol, 24.9, 20)
NHP_ENERGY_PRE[qol_ids] <- sim_nhp_dim(n_qol, 61.2, 25)
NHP_SOCIAL_PRE[qol_ids] <- sim_nhp_dim(n_qol, 18.4, 15)

# POST (better symptoms → lower NHP scores)
NHP_PAIN_POST[qol_ids]   <- sim_nhp_dim(n_qol, 14.5, 15)
NHP_EMO_POST[qol_ids]    <- sim_nhp_dim(n_qol, 8.1, 12)
NHP_MOB_POST[qol_ids]    <- sim_nhp_dim(n_qol, 13.3, 15)
NHP_SLEEP_POST[qol_ids]  <- sim_nhp_dim(n_qol, 26.5, 20)
NHP_ENERGY_POST[qol_ids] <- sim_nhp_dim(n_qol, 17.9, 20)
NHP_SOCIAL_POST[qol_ids] <- sim_nhp_dim(n_qol, 7.7, 10)

# EQ-5D utilities (here directly simulated around 0.30 pre, 0.60 post)
EQ5D_PRE  <- EQ5D_POST <- rep(NA_real_, n_patients)
EQ5D_PRE[qol_ids]  <- pmin(pmax(rnorm(n_qol, mean = 0.30, sd = 0.15), 0), 1)
EQ5D_POST[qol_ids] <- pmin(pmax(rnorm(n_qol, mean = 0.60, sd = 0.15), 0), 1)

############################################################
## 4. Assemble synthetic dataset ----
############################################################

df <- tibble(
  ID = ID,
  SEX = SEX,
  AGE = AGE,
  DIAGNOSIS = DIAGNOSIS,
  TRANSPLANTED = TRANSPLANTED,
  SECOND_TX = SECOND_TX,
  SUPER_URGENT = SUPER_URGENT,
  DELAI = DELAI,
  EVENT = EVENT,
  NHP_PAIN_PRE = NHP_PAIN_PRE,
  NHP_EMO_PRE = NHP_EMO_PRE,
  NHP_MOB_PRE = NHP_MOB_PRE,
  NHP_SLEEP_PRE = NHP_SLEEP_PRE,
  NHP_ENERGY_PRE = NHP_ENERGY_PRE,
  NHP_SOCIAL_PRE = NHP_SOCIAL_PRE,
  NHP_PAIN_POST = NHP_PAIN_POST,
  NHP_EMO_POST = NHP_EMO_POST,
  NHP_MOB_POST = NHP_MOB_POST,
  NHP_SLEEP_POST = NHP_SLEEP_POST,
  NHP_ENERGY_POST = NHP_ENERGY_POST,
  NHP_SOCIAL_POST = NHP_SOCIAL_POST,
  EQ5D_PRE = EQ5D_PRE,
  EQ5D_POST = EQ5D_POST
)

############################################################
## 5. Weibull survival modeling ----
############################################################

# Restrict to transplanted for post-TP survival
df_tp <- df %>% filter(TRANSPLANTED == 1)

fit_weib <- flexsurvreg(Surv(DELAI, EVENT) ~ 1, data = df_tp, dist = "weibull")
print(fit_weib)

shape_hat <- fit_weib$res["shape","est"]
scale_hat <- fit_weib$res["scale","est"]
cat("\nEstimated Weibull shape =", round(shape_hat, 3),
    "scale =", round(scale_hat, 1), "days\n")

# Survival at 1–5 years
S_weib <- function(t, k = shape_hat, la = scale_hat) exp(- (t / la)^k)
years <- 1:5
days  <- years * 365
S_tbl <- tibble(
  Year = years,
  Days = days,
  S    = S_weib(days)
)
print(S_tbl)

# Median and mean survival
median_surv <- flexsurv:::qgeneric(fit_weib, q = 0.5, type = "quantile")
mean_surv   <- scale_hat * gamma(1 + 1/shape_hat)

cat("\nMedian survival ≈", round(median_surv/365, 2), "years\n")
cat("Mean survival   ≈", round(mean_surv/365, 2), "years\n\n")

############################################################
## 6. Cost–Utility Analysis (base case) ----
############################################################

# Base HEOR assumptions
HR_nonTP_vs_TP <- 1.8        # worse survival without transplant
u_TP  <- 0.70
u_nonTP <- 0.60
c_TP_init <- 300000          # initial transplant cost (EUR)
c_TP_annual <- 20000         # annual post-TP cost
c_nonTP_annual <- 15000      # annual medical management cost
disc_rate <- 0.03
horizon_y <- 10
dt <- 30.4375                # monthly step (~days)

# Time grid and discounting
t_seq <- seq(0, horizon_y*365, by = dt)
disc  <- (1/(1+disc_rate))^(t_seq/365)

# Survival for TP and non-TP
S_TP    <- S_weib(t_seq)
S_nonTP <- S_TP^HR_nonTP_vs_TP

# QALYs
QALY_TP    <- sum(S_TP    * u_TP    * disc) * (dt/365)
QALY_nonTP <- sum(S_nonTP * u_nonTP * disc) * (dt/365)

# Costs
Cost_TP    <- c_TP_init + sum(S_TP    * (c_TP_annual/12)    * disc) * (dt/30.4375)
Cost_nonTP <- 0         + sum(S_nonTP * (c_nonTP_annual/12) * disc) * (dt/30.4375)

Delta_C <- Cost_TP - Cost_nonTP
Delta_E <- QALY_TP - QALY_nonTP
ICER    <- Delta_C / Delta_E

cat("Base-case QALY_TP    =", round(QALY_TP, 3), "\n")
cat("Base-case QALY_nonTP =", round(QALY_nonTP, 3), "\n")
cat("Base-case ΔCost      =", round(Delta_C, 0), "EUR\n")
cat("Base-case ΔQALY      =", round(Delta_E, 3), "\n")
cat("Base-case ICER       =", round(ICER, 0), "EUR/QALY\n\n")

############################################################
## 7. One-way sensitivity analysis (ICER) ----
############################################################

oneway <- tibble(
  param = c("HR", "u_TP", "u_nonTP", "c_TP_init", "c_TP_annual", "c_nonTP_annual"),
  base  = c(HR_nonTP_vs_TP, u_TP, u_nonTP, c_TP_init, c_TP_annual, c_nonTP_annual),
  low   = c(HR_nonTP_vs_TP*0.8, u_TP*0.9,     u_nonTP*0.9,
            c_TP_init*0.8,      c_TP_annual*0.8, c_nonTP_annual*0.8),
  high  = c(HR_nonTP_vs_TP*1.2, u_TP*1.1,     u_nonTP*1.1,
            c_TP_init*1.2,      c_TP_annual*1.2, c_nonTP_annual*1.2)
)

compute_ICER <- function(HR, u_tp, u_ng, c_init, c_tp, c_ng) {
  S_TP    <- S_weib(t_seq)
  S_nonTP <- S_TP^HR
  Q_tp <- sum(S_TP    * u_tp * disc) * (dt/365)
  Q_ng <- sum(S_nonTP * u_ng * disc) * (dt/365)
  C_tp <- c_init + sum(S_TP    * (c_tp/12) * disc) * (dt/30.4375)
  C_ng <- 0      + sum(S_nonTP * (c_ng/12) * disc) * (dt/30.4375)
  (C_tp - C_ng) / (Q_tp - Q_ng)
}

sens <- oneway %>%
  rowwise() %>%
  mutate(
    ICER_low = compute_ICER(
      HR     = if_else(param == "HR", low, HR_nonTP_vs_TP),
      u_tp   = if_else(param == "u_TP", low, u_TP),
      u_ng   = if_else(param == "u_nonTP", low, u_nonTP),
      c_init = if_else(param == "c_TP_init", low, c_TP_init),
      c_tp   = if_else(param == "c_TP_annual", low, c_TP_annual),
      c_ng   = if_else(param == "c_nonTP_annual", low, c_nonTP_annual)
    ),
    ICER_high = compute_ICER(
      HR     = if_else(param == "HR", high, HR_nonTP_vs_TP),
      u_tp   = if_else(param == "u_TP", high, u_TP),
      u_ng   = if_else(param == "u_nonTP", high, u_nonTP),
      c_init = if_else(param == "c_TP_init", high, c_TP_init),
      c_tp   = if_else(param == "c_TP_annual", high, c_TP_annual),
      c_ng   = if_else(param == "c_nonTP_annual", high, c_nonTP_annual)
    )
  ) %>%
  ungroup() %>%
  mutate(
    lo = pmin(ICER_low, ICER_high),
    hi = pmax(ICER_low, ICER_high)
  ) %>%
  arrange(desc(hi - lo))

print(sens)

############################################################
## 8. Probabilistic Sensitivity Analysis (PSA) ----
############################################################

nsim <- 1000

# Draw Weibull parameters (approx on log-scale)
vc  <- fit_weib$cov
par_mean <- c(log(scale_hat), log(shape_hat))
G <- diag(c(1/scale_hat, 1/shape_hat))
vc_log <- G %*% vc %*% t(G)
draws <- MASS::mvrnorm(nsim, mu = par_mean, Sigma = vc_log)
scale_s <- exp(draws[,1])
shape_s <- exp(draws[,2])

S_weib_s <- function(t, la, k) exp(- (t/la)^k)

# HR, utilities, costs
HR_s <- rlnorm(nsim, meanlog = log(HR_nonTP_vs_TP), sdlog = 0.2)

beta_mom <- function(mean, sd){
  v <- sd^2
  a <- ((1-mean)/v - 1/mean) * mean^2
  b <- a*(1/mean - 1)
  c(a = max(a, 1e-3), b = max(b, 1e-3))
}
ab_TP  <- beta_mom(u_TP, 0.05)
ab_NG  <- beta_mom(u_nonTP, 0.05)
u_TP_s  <- rbeta(nsim, ab_TP["a"], ab_TP["b"])
u_NG_s  <- rbeta(nsim, ab_NG["a"], ab_NG["b"])

gamma_from_mean_cv <- function(mean, cv){
  shape <- 1/(cv^2)
  scale <- mean * cv^2
  c(shape = shape, scale = scale)
}
g_init <- gamma_from_mean_cv(c_TP_init, 0.3)
g_tp   <- gamma_from_mean_cv(c_TP_annual, 0.3)
g_ng   <- gamma_from_mean_cv(c_nonTP_annual, 0.3)
c_init_s <- rgamma(nsim, shape = g_init["shape"], scale = g_init["scale"])
c_tp_s   <- rgamma(nsim, shape = g_tp["shape"],   scale = g_tp["scale"])
c_ng_s   <- rgamma(nsim, shape = g_ng["shape"],   scale = g_ng["scale"])

calc_ce <- function(la, k, HR, u_tp, u_ng, c_init, c_tp, c_ng){
  S_TP    <- S_weib_s(t_seq, la, k)
  S_nonTP <- S_TP^HR
  Q_tp <- sum(S_TP    * u_tp * disc) * (dt/365)
  Q_ng <- sum(S_nonTP * u_ng * disc) * (dt/365)
  C_tp <- c_init + sum(S_TP    * (c_tp/12) * disc) * (dt/30.4375)
  C_ng <- 0      + sum(S_nonTP * (c_ng/12) * disc) * (dt/30.4375)
  c(DeltaC = C_tp - C_ng, DeltaE = Q_tp - Q_ng)
}

psa_res <- purrr::pmap_dfr(
  list(scale_s, shape_s, HR_s, u_TP_s, u_NG_s, c_init_s, c_tp_s, c_ng_s),
  \(la, k, HR, u1, u0, ci, ct, cn){
    out <- calc_ce(la, k, HR, u1, u0, ci, ct, cn)
    tibble(DeltaC = out["DeltaC"], DeltaE = out["DeltaE"])
  }
) %>%
  mutate(ICER = DeltaC / DeltaE)

summary(psa_res$ICER)

cat("PSA completed: ", nrow(psa_res), "simulations\n")

############################################################
## End of script
############################################################
