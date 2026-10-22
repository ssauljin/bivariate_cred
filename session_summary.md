# Session summary — Poisson GLM vs. Poisson–gamma GLMM for regional loss counts

Handoff notes for continuing this analysis in a new session. Attach `loss_perregion_2015_2024.csv` and `poisson_models.R` alongside this file.

## User preferences
- **Code in R** (not Python). The current script uses base R only, with no packages.
- If a relevant skill is available, check it before answering.

## Data
- File: `loss_perregion_2015_2024.csv` (170 rows = 17 Korean regions × 10 years, 2015–2024; read with `fileEncoding = "UTF-8-BOM"`).
- Columns used: `Year`, `Region`, `AllCount` (response), `RecvBudget_m1Y` (recovery budget, lagged 1 year), `GRDPperbil`, `Population`, `Area`.
- Other columns that have not been used yet: AllLoss, ExtCount, ExtLoss, AvgAllLoss, PriorAllLoss, PriorExtLoss, AvgExtLoss, RecvBudget, RecvBudget_m2Y, RecvBudget_m3Y, RecvBudget_mvAvg3Y, Freq_sanity, Loss_sanity.
- Training set: 2015–2023 (153 rows, 9 years per region). Test set: 2024 (17 rows).
- AllCount on the training set: mean 6.09, variance 18.2 (overdispersed).

## Models
- **Poisson GLM**: `glm(AllCount ~ log(covariates), family = poisson)`.
- **Poisson–gamma GLMM**: yᵢⱼ | uᵢ ~ Poisson(uᵢ μᵢⱼ), log μᵢⱼ = xᵢⱼ'β, uᵢ ~ Gamma(φ, φ), so E[uᵢ] = 1 and Var[uᵢ] = 1/φ (random effect on Region).
  - Fitted by **exact marginal ML** using the closed-form region likelihood:
    Lᵢ = Πⱼ μᵢⱼ^yᵢⱼ / yᵢⱼ! · φ^φ/Γ(φ) · Γ(φ+Yᵢ)/(φ+Mᵢ)^(φ+Yᵢ), where Yᵢ = Σⱼ yᵢⱼ and Mᵢ = Σⱼ μᵢⱼ.
  - Optimisation: BFGS (`optim`) with an analytic gradient, then a Newton polish. Parameterised as (β, log φ).
  - Standard errors come from the observed information (central-difference Hessian of the analytic gradient). SE of Var(u) uses the delta method.
  - Check: `integrate()` over each region reproduces the closed-form log-likelihood to 1e-6.
  - LR test for the random effect uses the boundary-corrected ½χ²₁ distribution.
- **Prediction for 2024 (T+1)**:
  - GLM: ŷ = μᵢ,T+1.
  - GLMM: ŷ = μᵢ,T+1 · E[uᵢ | yᵢ,1..T], where E[uᵢ | y] = (φ̂ + Yᵢ)/(φ̂ + M̂ᵢ), computed from training years only.
  - Criteria: RMSE, Poisson deviance 2Σ[y log(y/ŷ) − (y − ŷ)] (with 0·log 0 = 0), plus MAE for reference.

## Results

### Spec A — full: log(RecvBudget_m1Y) + log(GRDPperbil) + log(Population) + log(Area)

| Term | GLM est (SE) | GLMM est (SE) |
|---|---:|---:|
| Intercept | −3.192 (1.024) | −0.096 (2.492) |
| log(RecvBudget_m1Y) | 0.060 (0.067) | 0.033 (0.111) |
| log(GRDPperbil) | 0.313 (0.130) | 1.120 (0.432) |
| log(Population) | −0.117 (0.134) | −0.910 (0.415) |
| log(Area) | 0.247 (0.047) | 0.230 (0.088) |
| Var(uᵢ) | — | 0.091 (0.066); φ̂ = 10.98 |
| logLik / AIC / BIC | −388.78 / 787.55 / 802.71 | −382.69 / 777.38 / 795.56 |

- GLM Pearson χ²/df = 1.645.
- LR test for the random effect = 12.18 (p = 0.00024).
- 2024 prediction (GLM vs GLMM):
  - RMSE: 4.869 vs **4.561**
  - Deviance: 46.53 vs **41.01**
  - MAE: 3.49 vs 3.40
  - Predicted total: 113.0 vs 133.4 (actual 151)

### Spec B — no budget: log(GRDPperbil) + log(Population) + log(Area)

| Term | GLM est (SE) | GLMM est (SE) |
|---|---:|---:|
| Intercept | −2.559 (0.737) | 0.399 (1.939) |
| log(GRDPperbil) | 0.327 (0.129) | 1.179 (0.390) |
| log(Population) | −0.115 (0.133) | −0.958 (0.389) |
| log(Area) | 0.283 (0.026) | 0.247 (0.068) |
| Var(uᵢ) | — | 0.100 (0.064); φ̂ = 10.05 |
| logLik / AIC / BIC | −389.19 / 786.38 / 798.50 | −382.73 / 775.47 / 790.62 |

- LR test for the random effect = 12.91 (p = 0.00016).
- LR test for dropping the budget term: 0.82 (p = 0.37) in the GLM and 0.09 (p = 0.77) in the GLMM.
- 2024 prediction (GLM vs GLMM):
  - RMSE: 4.891 vs **4.554**
  - Deviance: 46.46 vs **40.87**
  - MAE: 3.49 vs 3.40
  - Predicted total: 111.6 vs 133.9 (actual 151)

## Key interpretations so far
- The region random effect is clearly supported: the GLMM wins on AIC/BIC and the LR test, and it predicts 2024 better on both criteria in both specs.
- The lagged recovery budget has no detectable effect. Dropping it improves AIC/BIC and leaves prediction essentially unchanged, which makes Spec B the preferred specification.
- log(Area) is the most stable effect across models (≈ +0.23 to +0.28).
- The GRDP and Population coefficients are unstable: they shift from 0.31/−0.12 in the GLM to about 1.1–1.2/−0.9 in the GLMM. Pooled log(GRDP) and log(Population) are correlated at r = 0.95, and Population barely varies within a region (within-region SD of the log ≈ 0.05). This is collinearity, and the pair effectively acts as a GRDP-per-capita effect.
- GLM model-based SEs understate uncertainty. Cluster-robust SEs are larger, for example log(Area) 0.074 vs 0.047.
- Empirical-Bayes multipliers range from Ulsan ≈ 0.53–0.56 and Seoul ≈ 0.61 up to Incheon ≈ 1.43–1.45.
- Both models under-predict the 2024 total. The largest misses are Gangwon (22 actual vs ≈ 9–12 predicted), Gyeongnam (19 vs ≈ 10) and Jeonnam (17 vs ≈ 9). This points to a year-level shock that the model has no term for (no year or trend effect).

## Code
- `poisson_models.R` defines `run_spec(covars, tag, train_years = 2015:2023, test_year = 2024)`.
  - It prints the regression table, fit statistics, empirical-Bayes effects, 2024 predictions and criteria.
  - It writes `regression_table_<tag>.csv`, `goodness_of_fit_<tag>.csv`, `region_random_effects_<tag>.csv`, `predictions_2024_<tag>.csv` and `prediction_criteria_2024_<tag>.csv`.
  - The bottom of the script currently runs `run_spec(c("RecvBudget_m1Y","GRDPperbil","Population","Area"), "full")` and `run_spec(c("GRDPperbil","Population","Area"), "nobudget")`.
  - A new specification is a one-line call. Note that the helper `design()` applies `log()` to every covariate, so a ratio such as GRDP per capita must first be added as a column in `df`.
- `poisson_models.py` is an earlier Python version of Spec A only. It is superseded by the R script.

## Suggested next steps (not yet done)
1. Replace log(GRDPperbil) and log(Population) with log(GRDPperbil/Population) to resolve the collinearity.
2. Add a linear year term or a post-2020 indicator to capture common time shocks and the 2024 under-prediction.
3. Fit a log-normal random-effect GLMM (`lme4::glmer`, nAGQ = 25) for comparison, and/or cross-check the gamma model with `hglm`.
4. Consider rolling-origin validation (for example, predicting 2022, 2023 and 2024 in turn), because 17 test points is a small sample.
