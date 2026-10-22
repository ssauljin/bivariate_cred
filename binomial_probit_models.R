## =============================================================================
## Extreme claim counts M_it given total claim counts N_it:
## binomial-probit GLM vs. binomial-probit GLMM (Gaussian random effect on Region).
## Estimation on 2015-2023, out-of-sample prediction of 2024.
##
## GLMM:  M_ij | N_ij, v_i ~ Bin(N_ij, Phi(eta_ij + v_i)),  eta_ij = x_ij' beta,
##        v_i ~ N(0, sigma^2).
## The region-level marginal likelihood is a one-dimensional integral,
##   L_i = prod_j choose(N_ij, M_ij)
##         * int phi(z) prod_j Phi(eta_ij + sigma z)^M_ij Phi(-eta_ij - sigma z)^(N_ij - M_ij) dz,
## evaluated by Gauss-Hermite quadrature (checked against integrate()).
## Periods with N_ij = 0 contribute nothing to the likelihood.
##
## The total claim count part N_ij | u_i ~ Pois(u_i mu_ij), u_i ~ Gamma(phi, phi), is
## fitted as in poisson_models.R; by the likelihood factorization of the thinning
## model the two parts are estimated separately.
##
## Prediction for T+1 = 2024:
##   p_hat_i = E[Phi(eta_i,T+1 + v_i) | history]                (posterior predictive)
##   (a) given N_i,T+1 :  M_hat = N_i,T+1 * p_hat_i
##   (b) unconditional :  M_hat = E[N_i,T+1 | history] * p_hat_i
##   GLM benchmark     :  p_hat_i = Phi(eta_i,T+1), E[N_i,T+1] = mu_i,T+1 (Poisson GLM)
##
## Base R only -- no packages required.  Run:  Rscript binomial_probit_models.R
## =============================================================================

df <- read.csv("loss_perregion_2015_2024.csv", fileEncoding = "UTF-8-BOM")
out_dir <- "results_extreme"
dir.create(out_dir, showWarnings = FALSE)

## ---- Helpers -----------------------------------------------------------------
design <- function(dat, covars) {
  X <- cbind(1, sapply(covars, function(v) log(dat[[v]])))
  colnames(X) <- c("(Intercept)", paste0("log(", covars, ")"))
  X
}

## Gauss-Hermite rule for N(0,1) (Golub-Welsch): sum_q w_q f(z_q) ~ E[f(Z)]
gh_normal <- function(Q) {
  J <- matrix(0, Q, Q)
  J[cbind(1:(Q - 1), 2:Q)] <- J[cbind(2:Q, 1:(Q - 1))] <- sqrt(1:(Q - 1))
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}

num_hessian <- function(f_grad, theta, h = 1e-5) {          # central diff of analytic gradient
  k <- length(theta); H <- matrix(0, k, k)
  for (j in seq_len(k)) { e <- rep(0, k); e[j] <- h
    H[, j] <- (f_grad(theta + e) - f_grad(theta - e)) / (2 * h) }
  (H + t(H)) / 2
}

newton_polish <- function(theta, grad_f) {
  for (it in 1:20) {
    g <- grad_f(theta); if (max(abs(g)) < 1e-9) break
    theta <- theta - solve(num_hessian(grad_f, theta), g)
  }
  theta
}

## ---- Poisson-gamma GLMM for N (same likelihood as poisson_models.R) ----------
fit_poisson_gamma <- function(y, X, grp, start) {
  p <- ncol(X); lgam_y <- lgamma(y + 1); Y_g <- as.numeric(tapply(y, grp, sum))
  nll <- function(theta) {
    beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
    eta <- as.vector(X %*% beta); M_g <- as.numeric(tapply(exp(eta), grp, sum))
    -(sum(y * eta - lgam_y) +
      sum(phi * logphi - lgamma(phi) + lgamma(phi + Y_g) - (phi + Y_g) * log(phi + M_g)))
  }
  gr <- function(theta) {
    beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
    mu <- as.vector(exp(X %*% beta)); M_g <- as.numeric(tapply(mu, grp, sum))
    w_g <- (phi + Y_g) / (phi + M_g)
    g_beta <- crossprod(X, y) - crossprod(X, mu * w_g[grp])
    g_lphi <- phi * sum(logphi + 1 - digamma(phi) + digamma(phi + Y_g) - log(phi + M_g) - w_g)
    -c(as.vector(g_beta), g_lphi)
  }
  opt <- optim(start, nll, gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-14))
  theta <- newton_polish(opt$par, gr)
  beta <- theta[1:p]; phi <- exp(theta[p + 1])
  M_g <- as.numeric(tapply(as.vector(exp(X %*% beta)), grp, sum))
  list(beta = beta, phi = phi, loglik = -nll(theta), u_post = (phi + Y_g) / (phi + M_g))
}

## ---- Main routine --------------------------------------------------------------
run_spec <- function(covars, tag, train_years = 2015:2023, test_year = 2024, Q = 100) {

  cat("\n\n#################################################################\n")
  cat("## Specification:", tag, "\n## Covariates :", paste0("log(", covars, ")", collapse = " + "), "\n")
  cat("#################################################################\n")

  d <- subset(df, Year %in% train_years); d <- d[order(d$Region, d$Year), ]
  d$Region <- factor(d$Region)
  X_all <- design(d, covars); grp_all <- as.integer(d$Region); G <- nlevels(d$Region)

  ## binomial part uses the region-years with at least one claim
  keep <- d$AllCount > 0
  n <- d$AllCount[keep]; m <- d$ExtCount[keep]; X <- X_all[keep, , drop = FALSE]
  grp <- grp_all[keep]; p <- ncol(X); nobs <- sum(keep)
  stopifnot(length(unique(grp)) == G)   # every region has at least one claim
  stopifnot(all(m <= n))
  cat(sprintf("Binomial part: %d region-years with N > 0 (of %d), %d claims, %d extreme\n",
              nobs, nrow(d), sum(n), sum(m)))

  ## 1. Probit GLM ---------------------------------------------------------------
  db <- d[keep, ]
  fml <- as.formula(paste("cbind(ExtCount, AllCount - ExtCount) ~",
                          paste0("log(", covars, ")", collapse = " + ")))
  glm_fit <- glm(fml, family = binomial(link = "probit"), data = db)
  glm_ll  <- as.numeric(logLik(glm_fit)); glm_k <- p
  glm_aic <- -2 * glm_ll + 2 * glm_k; glm_bic <- -2 * glm_ll + log(nobs) * glm_k
  pearson_disp <- sum(residuals(glm_fit, type = "pearson")^2) / glm_fit$df.residual

  ## 2. Probit GLMM, marginal ML by Gauss-Hermite; theta = (beta, log sigma) -----
  gh <- gh_normal(Q); lchoose_nm <- sum(lchoose(n, m))

  ## Q x nobs matrix of per-observation conditional log-likelihoods at each node
  node_terms <- function(beta, sigma) {
    eta <- as.vector(X %*% beta)
    A <- outer(sigma * gh$z, eta, "+")                         # a_qj = eta_j + sigma z_q
    list(A = A, ll = sweep(pnorm(A, log.p = TRUE), 2, m, "*") +
                     sweep(pnorm(-A, log.p = TRUE), 2, n - m, "*"))
  }
  ## Q x G matrix: log w_q + sum_{j in region} ll_qj
  region_logpost <- function(ll) {
    S <- t(rowsum(t(ll), grp, reorder = TRUE))
    sweep(S, 1, log(gh$w), "+")
  }
  logsumexp_cols <- function(L) { mx <- apply(L, 2, max); mx + log(colSums(exp(sweep(L, 2, mx)))) }

  nll <- function(theta) {
    nt <- node_terms(theta[1:p], exp(theta[p + 1]))
    -(lchoose_nm + sum(logsumexp_cols(region_logpost(nt$ll))))
  }
  gr <- function(theta) {
    beta <- theta[1:p]; sigma <- exp(theta[p + 1])
    nt <- node_terms(beta, sigma); L <- region_logpost(nt$ll)
    P <- exp(sweep(L, 2, logsumexp_cols(L)))                    # Q x G posterior node weights
    ## d ll_qj / d a_qj = m_j phi/Phi(a) - (n_j - m_j) phi/Phi(-a)
    mills <- function(a) exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
    D <- sweep(mills(nt$A), 2, m, "*") - sweep(mills(-nt$A), 2, n - m, "*")
    Pobs <- P[, grp, drop = FALSE]                               # Q x nobs
    g_beta <- crossprod(X, colSums(Pobs * D))
    g_lsig <- sigma * sum(colSums(Pobs * D * gh$z))
    -c(as.vector(g_beta), g_lsig)
  }

  opt <- optim(c(coef(glm_fit), log(0.5)), nll, gr, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-14))
  theta_hat <- newton_polish(opt$par, gr)
  cat(sprintf("GLMM converged | max|grad| = %.2e\n", max(abs(gr(theta_hat)))))

  V  <- solve(num_hessian(gr, theta_hat)); se <- sqrt(diag(V))
  names(theta_hat) <- names(se) <- c(colnames(X), "log(sigma)")
  beta_hat <- theta_hat[1:p]; sigma_hat <- exp(theta_hat[p + 1])
  glmm_ll  <- -nll(theta_hat); glmm_k <- p + 1
  glmm_aic <- -2 * glmm_ll + 2 * glmm_k; glmm_bic <- -2 * glmm_ll + log(nobs) * glmm_k
  var_v    <- sigma_hat^2; se_var_v <- 2 * var_v * se[p + 1]     # delta method
  ci_var_v <- exp(2 * (theta_hat[p + 1] + c(-1.96, 1.96) * se[p + 1]))

  ## quadrature checks: more nodes, and adaptive integrate() per region
  ll_200 <- { gh_save <- gh; gh <- gh_normal(200); v <- -nll(theta_hat); gh <- gh_save; v }
  eta_hat <- as.vector(X %*% beta_hat)
  num_ll <- lchoose_nm + sum(sapply(seq_len(G), function(g) {
    idx <- grp == g
    if (!any(idx)) return(0)
    f <- function(z) sapply(z, function(zz) exp(sum(m[idx] * pnorm(eta_hat[idx] + sigma_hat * zz, log.p = TRUE) +
                                                  (n[idx] - m[idx]) * pnorm(-eta_hat[idx] - sigma_hat * zz, log.p = TRUE)))) * dnorm(z)
    log(integrate(f, -10, 10, rel.tol = 1e-12, subdivisions = 2000L)$value)
  }))
  cat(sprintf("loglik: GH(%d) = %.6f | GH(200) = %.6f | integrate() = %.6f\n",
              Q, glmm_ll, ll_200, num_ll))

  lr <- 2 * (glmm_ll - glm_ll); lr_p <- 0.5 * pchisq(lr, 1, lower.tail = FALSE)

  ## regression table -----------------------------------------------------------
  glm_se <- sqrt(diag(vcov(glm_fit))); glm_z <- coef(glm_fit) / glm_se; glmm_z <- beta_hat / se[1:p]
  tab <- data.frame(Term = colnames(X),
                    GLM_Estimate = coef(glm_fit), GLM_SE = glm_se, GLM_z = glm_z,
                    GLM_p = 2 * pnorm(-abs(glm_z)),
                    GLMM_Estimate = beta_hat, GLMM_SE = se[1:p], GLMM_z = glmm_z,
                    GLMM_p = 2 * pnorm(-abs(glmm_z)), row.names = NULL)
  tab <- rbind(tab, data.frame(Term = "Var(v_i) = sigma^2 [Region]",
                               GLM_Estimate = NA, GLM_SE = NA, GLM_z = NA, GLM_p = NA,
                               GLMM_Estimate = var_v, GLMM_SE = se_var_v, GLMM_z = NA, GLMM_p = NA))
  fit <- data.frame(
    Statistic = c("N (region-years with N>0)", "Regions", "Parameters", "Log-likelihood",
                  "AIC", "BIC", "Pearson chi2 / df"),
    Probit_GLM  = c(nobs, G, glm_k, glm_ll, glm_aic, glm_bic, pearson_disp),
    Probit_GLMM = c(nobs, G, glmm_k, glmm_ll, glmm_aic, glmm_bic, NA))

  cat("\n=== Coefficients ===\n");    print(tab, digits = 4, row.names = FALSE)
  cat("\n=== Goodness of fit ===\n"); print(fit, digits = 6, row.names = FALSE)
  cat(sprintf("\nsigma_hat = %.4f, Var(v) = %.4f (95%% CI %.4f - %.4f)\n",
              sigma_hat, var_v, ci_var_v[1], ci_var_v[2]))
  cat(sprintf("LR test for random effect: LR = %.3f, p (0.5*chi2_1) = %.3g\n", lr, lr_p))

  ## posterior of v_i (Gauss-Hermite with Q = 200) --------------------------------
  gh <- gh_normal(200)
  nt <- node_terms(beta_hat, sigma_hat); L <- region_logpost(nt$ll)
  P  <- exp(sweep(L, 2, logsumexp_cols(L)))                       # Q x G
  vq <- sigma_hat * gh$z
  v_post <- as.vector(crossprod(P, vq))
  ## regions with no claims in the training years keep the prior (all-zero columns of ll)
  Y_g <- as.numeric(tapply(d$AllCount, grp_all, sum)); M_g <- as.numeric(tapply(d$ExtCount, grp_all, sum))

  ## 3. Out-of-sample prediction for test_year ----------------------------------
  d24 <- subset(df, Year == test_year); d24 <- d24[order(d24$Region), ]
  stopifnot(identical(as.character(d24$Region), levels(d$Region)))
  X24 <- design(d24, covars); N24 <- d24$AllCount; M24 <- d24$ExtCount
  eta24 <- as.vector(X24 %*% beta_hat)
  p_hat_glmm <- as.vector(colSums(P * pnorm(outer(vq, eta24, "+"))))   # E[Phi(eta + v) | H]
  p_prior    <- pnorm(eta24 / sqrt(1 + var_v))
  p_hat_glm  <- as.vector(pnorm(X24 %*% coef(glm_fit)))

  ## total claim count part, same covariates (Poisson GLM and Poisson-gamma GLMM)
  yN <- d$AllCount
  pois_glm <- glm(as.formula(paste("AllCount ~", paste0("log(", covars, ")", collapse = " + "))),
                  family = poisson(link = "log"), data = d)
  pg <- fit_poisson_gamma(yN, X_all, grp_all, c(coef(pois_glm), log(2)))
  EN_glm  <- as.vector(exp(X24 %*% coef(pois_glm)))
  EN_glmm <- as.vector(exp(X24 %*% pg$beta)) * pg$u_post
  cat(sprintf("\nN part (Poisson-gamma GLMM): phi_hat = %.4f, loglik = %.4f\n", pg$phi, pg$loglik))

  pred_tab <- data.frame(Region = d24$Region, N_train = Y_g, M_train = M_g,
                         E_v_given_H = v_post, p_prior = p_prior,
                         p_hat_GLM = p_hat_glm, p_hat_GLMM = p_hat_glmm,
                         N_2024 = N24, M_2024 = M24,
                         M_given_N_GLM = N24 * p_hat_glm, M_given_N_GLMM = N24 * p_hat_glmm,
                         EN_GLM = EN_glm, EN_GLMM = EN_glmm,
                         M_uncond_GLM = EN_glm * p_hat_glm, M_uncond_GLMM = EN_glmm * p_hat_glmm)

  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
  mae  <- function(y, yhat) mean(abs(y - yhat))
  xlogx_ratio <- function(a, b) ifelse(a > 0, a * log(a / b), 0)
  binom_dev <- function(m, n, p) 2 * sum(xlogx_ratio(m, n * p) + xlogx_ratio(n - m, n * (1 - p)))
  pois_dev  <- function(y, yhat) 2 * sum(xlogx_ratio(y, yhat) - (y - yhat))
  crit <- data.frame(
    Criterion = c("(a) given N: RMSE", "(a) given N: MAE", "(a) given N: binomial deviance",
                  "(b) unconditional: RMSE", "(b) unconditional: MAE", "(b) unconditional: Poisson deviance",
                  "Predicted total, given N", "Predicted total, unconditional"),
    GLM  = c(rmse(M24, N24 * p_hat_glm), mae(M24, N24 * p_hat_glm), binom_dev(M24, N24, p_hat_glm),
             rmse(M24, EN_glm * p_hat_glm), mae(M24, EN_glm * p_hat_glm), pois_dev(M24, EN_glm * p_hat_glm),
             sum(N24 * p_hat_glm), sum(EN_glm * p_hat_glm)),
    GLMM = c(rmse(M24, N24 * p_hat_glmm), mae(M24, N24 * p_hat_glmm), binom_dev(M24, N24, p_hat_glmm),
             rmse(M24, EN_glmm * p_hat_glmm), mae(M24, EN_glmm * p_hat_glmm), pois_dev(M24, EN_glmm * p_hat_glmm),
             sum(N24 * p_hat_glmm), sum(EN_glmm * p_hat_glmm)))

  cat(sprintf("\n=== %d hold-out predictions by region ===\n", test_year))
  print(pred_tab[order(-pred_tab$E_v_given_H), ], digits = 3, row.names = FALSE)
  cat(sprintf("\n=== %d prediction criteria (n = %d regions); actual total M = %d ===\n",
              test_year, G, sum(M24)))
  print(crit, digits = 4, row.names = FALSE)

  ## write outputs -------------------------------------------------------------------
  write.csv(tab,      file.path(out_dir, paste0("regression_table_M_", tag, ".csv")),   row.names = FALSE)
  write.csv(fit,      file.path(out_dir, paste0("goodness_of_fit_M_", tag, ".csv")),    row.names = FALSE)
  write.csv(pred_tab, file.path(out_dir, paste0("predictions_M_", test_year, "_", tag, ".csv")), row.names = FALSE)
  write.csv(crit,     file.path(out_dir, paste0("prediction_criteria_M_", test_year, "_", tag, ".csv")), row.names = FALSE)
  invisible(list(table = tab, fit = fit, pred = pred_tab, crit = crit, lr = lr, lr_p = lr_p,
                 sigma = sigma_hat, glm = glm_fit))
}

## ---- Descriptive: share of extreme claims by year -------------------------------
by_year <- aggregate(cbind(AllCount, ExtCount) ~ Year, data = df, FUN = sum)
by_year$Share <- by_year$ExtCount / by_year$AllCount
cat("=== Extreme claims by year ===\n"); print(by_year, digits = 3, row.names = FALSE)
write.csv(by_year, file.path(out_dir, "extreme_share_by_year.csv"), row.names = FALSE)

## ---- Run both specifications ----------------------------------------------------
full     <- run_spec(c("RecvBudget_m1Y", "GRDPperbil", "Population", "Area"), tag = "full")
nobudget <- run_spec(c("GRDPperbil", "Population", "Area"),                   tag = "nobudget")

## ---- Diagnostic: year and region fixed effects in the probit GLM (Spec B) ----------
dfe <- subset(df, Year <= 2023 & AllCount > 0)
fe0 <- glm(cbind(ExtCount, AllCount - ExtCount) ~ log(GRDPperbil) + log(Population) + log(Area),
           family = binomial(link = "probit"), data = dfe)
fe1 <- update(fe0, . ~ . + factor(Year))
fe2 <- update(fe1, . ~ . + factor(Region))
pdisp <- function(f) sum(residuals(f, type = "pearson")^2) / f$df.residual
cat("\n=== Fixed-effect diagnostic (probit GLM, Spec B covariates) ===\n")
cat(sprintf("Pearson chi2/df: base %.2f | + year FE %.2f | + year and region FE %.2f\n",
            pdisp(fe0), pdisp(fe1), pdisp(fe2)))
print(anova(fe0, fe1, fe2, test = "Chisq"))
