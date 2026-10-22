## =============================================================================
## Poisson GLM vs. Poisson-gamma GLMM (gamma random effect on Region, mean 1)
## Estimation on 2015-2023, out-of-sample prediction of 2024.
##
## GLMM:  y_ij | u_i ~ Poisson(u_i * mu_ij),  log(mu_ij) = x_ij' beta,
##        u_i ~ Gamma(shape = phi, rate = phi)  =>  E[u_i] = 1, Var[u_i] = 1/phi.
## Integrating u_i out gives a closed-form region-level marginal likelihood
##   L_i = prod_j mu_ij^y_ij / y_ij!  *  phi^phi / Gamma(phi)
##         * Gamma(phi + Y_i) / (phi + M_i)^(phi + Y_i),   Y_i = sum_j y_ij, M_i = sum_j mu_ij,
## so the GLMM is fitted by EXACT maximum likelihood (no Laplace / quadrature).
##
## Prediction for T+1 = 2024:
##   GLM : yhat_i = mu_i,T+1 = exp(x_i,T+1' beta_GLM)
##   GLMM: yhat_i = mu_i,T+1 * E[u_i | y_i1..y_iT],  E[u_i | y] = (phi + Y_i)/(phi + M_i)
##
## Base R only -- no packages required.  Run:  Rscript poisson_models.R
## =============================================================================

df <- read.csv("loss_perregion_2015_2024.csv", fileEncoding = "UTF-8-BOM")

## ---- Helper: build a design matrix from a vector of covariate names ---------
design <- function(dat, covars) {
  X <- cbind(1, sapply(covars, function(v) log(dat[[v]])))
  colnames(X) <- c("(Intercept)", paste0("log(", covars, ")"))
  X
}

## ---- Main routine: fit both models on 2015-2023, predict 2024 --------------
run_spec <- function(covars, tag, train_years = 2015:2023, test_year = 2024) {

  cat("\n\n#################################################################\n")
  cat("## Specification:", tag, "\n## Covariates :", paste0("log(", covars, ")", collapse = " + "), "\n")
  cat("#################################################################\n")

  d <- subset(df, Year %in% train_years); d <- d[order(d$Region, d$Year), ]
  d$Region <- factor(d$Region)
  y   <- d$AllCount
  X   <- design(d, covars)
  grp <- as.integer(d$Region); G <- nlevels(d$Region); n <- length(y); p <- ncol(X)
  cat(sprintf("n = %d obs, %d regions, years %d-%d\n", n, G, min(d$Year), max(d$Year)))

  ## 1. Poisson GLM -----------------------------------------------------------
  fml <- as.formula(paste("AllCount ~", paste0("log(", covars, ")", collapse = " + ")))
  glm_fit <- glm(fml, family = poisson(link = "log"), data = d)
  glm_ll  <- as.numeric(logLik(glm_fit)); glm_k <- p
  glm_aic <- -2 * glm_ll + 2 * glm_k;     glm_bic <- -2 * glm_ll + log(n) * glm_k
  pearson_disp  <- sum(residuals(glm_fit, type = "pearson")^2) / glm_fit$df.residual
  deviance_disp <- deviance(glm_fit) / glm_fit$df.residual

  ## 2. Poisson-gamma GLMM, exact marginal ML; theta = (beta, log phi) --------
  lgam_y <- lgamma(y + 1)
  Y_g    <- as.numeric(tapply(y, grp, sum))

  negloglik <- function(theta) {
    beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
    eta  <- as.vector(X %*% beta); mu <- exp(eta)
    M_g  <- as.numeric(tapply(mu, grp, sum))
    -(sum(y * eta - lgam_y) +
      sum(phi * logphi - lgamma(phi) + lgamma(phi + Y_g) - (phi + Y_g) * log(phi + M_g)))
  }
  grad_negloglik <- function(theta) {
    beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
    mu   <- as.vector(exp(X %*% beta))
    M_g  <- as.numeric(tapply(mu, grp, sum))
    w_g  <- (phi + Y_g) / (phi + M_g)                         # E[u_i | y]
    g_beta <- crossprod(X, y) - crossprod(X, mu * w_g[grp])
    g_lphi <- phi * sum(logphi + 1 - digamma(phi) + digamma(phi + Y_g) - log(phi + M_g) - w_g)
    -c(as.vector(g_beta), g_lphi)
  }
  num_hessian <- function(f_grad, theta, h = 1e-5) {         # central diff of analytic gradient
    k <- length(theta); H <- matrix(0, k, k)
    for (j in seq_len(k)) { e <- rep(0, k); e[j] <- h
      H[, j] <- (f_grad(theta + e) - f_grad(theta - e)) / (2 * h) }
    (H + t(H)) / 2
  }

  opt <- optim(c(coef(glm_fit), log(2)), negloglik, grad_negloglik, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-14))
  theta_hat <- opt$par
  for (it in 1:20) {                                          # Newton polish
    g <- grad_negloglik(theta_hat); if (max(abs(g)) < 1e-9) break
    theta_hat <- theta_hat - solve(num_hessian(grad_negloglik, theta_hat), g)
  }
  cat(sprintf("GLMM converged | max|grad| = %.2e\n", max(abs(grad_negloglik(theta_hat)))))

  V  <- solve(num_hessian(grad_negloglik, theta_hat))         # inverse observed information
  se <- sqrt(diag(V))
  names(theta_hat) <- names(se) <- c(colnames(X), "log(phi)")

  glmm_ll  <- -negloglik(theta_hat); glmm_k <- p + 1
  glmm_aic <- -2 * glmm_ll + 2 * glmm_k; glmm_bic <- -2 * glmm_ll + log(n) * glmm_k
  beta_hat <- theta_hat[1:p]; phi_hat <- exp(theta_hat[p + 1])
  var_u    <- 1 / phi_hat; se_var_u <- var_u * se[p + 1]     # delta method
  ci_var_u <- exp(-(theta_hat[p + 1] + c(1.96, -1.96) * se[p + 1]))

  ## quadrature check of the closed-form marginal likelihood
  num_ll <- sum(sapply(seq_len(G), function(g) {
    idx <- grp == g; mu <- exp(X[idx, , drop = FALSE] %*% beta_hat); yg <- y[idx]
    logf <- function(uu) sum(dpois(yg, uu * mu, log = TRUE)) + dgamma(uu, phi_hat, phi_hat, log = TRUE)
    m <- logf(1); f <- function(u) sapply(u, function(uu) exp(logf(uu) - m))
    a <- phi_hat + sum(yg); b <- phi_hat + sum(mu)
    m + log(integrate(f, qgamma(1e-14, a, b), qgamma(1 - 1e-14, a, b),
                      rel.tol = 1e-12, subdivisions = 1000L)$value)
  }))
  cat(sprintf("closed-form loglik = %.6f | numerical quadrature = %.6f\n", glmm_ll, num_ll))

  lr <- 2 * (glmm_ll - glm_ll); lr_p <- 0.5 * pchisq(lr, 1, lower.tail = FALSE)

  ## regression table ---------------------------------------------------------
  glm_se <- sqrt(diag(vcov(glm_fit))); glm_z <- coef(glm_fit) / glm_se; glmm_z <- beta_hat / se[1:p]
  tab <- data.frame(Term = colnames(X),
                    GLM_Estimate = coef(glm_fit), GLM_SE = glm_se, GLM_z = glm_z,
                    GLM_p = 2 * pnorm(-abs(glm_z)),
                    GLMM_Estimate = beta_hat, GLMM_SE = se[1:p], GLMM_z = glmm_z,
                    GLMM_p = 2 * pnorm(-abs(glmm_z)), row.names = NULL)
  tab <- rbind(tab, data.frame(Term = "Var(u_i) = 1/phi [Region]",
                               GLM_Estimate = NA, GLM_SE = NA, GLM_z = NA, GLM_p = NA,
                               GLMM_Estimate = var_u, GLMM_SE = se_var_u, GLMM_z = NA, GLMM_p = NA))
  fit <- data.frame(
    Statistic = c("N (obs)", "Regions", "Parameters", "Log-likelihood", "AIC", "BIC",
                  "Deviance", "Pearson chi2 / df", "Deviance / df"),
    Poisson_GLM        = c(n, G, glm_k, glm_ll, glm_aic, glm_bic, deviance(glm_fit), pearson_disp, deviance_disp),
    Poisson_gamma_GLMM = c(n, G, glmm_k, glmm_ll, glmm_aic, glmm_bic, NA, NA, NA))

  cat("\n=== Coefficients ===\n");    print(tab, digits = 4, row.names = FALSE)
  cat("\n=== Goodness of fit ===\n"); print(fit, digits = 6, row.names = FALSE)
  cat(sprintf("\nphi_hat = %.4f, Var(u) = %.4f (95%% CI %.4f - %.4f)\n", phi_hat, var_u, ci_var_u[1], ci_var_u[2]))
  cat(sprintf("LR test for random effect: LR = %.3f, p (0.5*chi2_1) = %.3g\n", lr, lr_p))

  ## empirical-Bayes region effects -------------------------------------------
  mu_hat <- as.vector(exp(X %*% beta_hat)); M_g <- as.numeric(tapply(mu_hat, grp, sum))
  u_post <- (phi_hat + Y_g) / (phi_hat + M_g)
  eb <- data.frame(Region = levels(d$Region), sum_y = Y_g, sum_mu = M_g, EB_u_i = u_post)
  cat("\n=== Empirical Bayes region effects u_i ===\n")
  print(eb[order(-eb$EB_u_i), ], digits = 4, row.names = FALSE)

  ## 3. Out-of-sample prediction for test_year ---------------------------------
  d24 <- subset(df, Year == test_year); d24 <- d24[order(d24$Region), ]
  stopifnot(identical(as.character(d24$Region), levels(d$Region)))
  X24 <- design(d24, covars); y24 <- d24$AllCount
  pred_glm  <- as.vector(exp(X24 %*% coef(glm_fit)))
  mu24_glmm <- as.vector(exp(X24 %*% beta_hat))
  pred_glmm <- mu24_glmm * u_post

  rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
  poisson_deviance <- function(y, yhat) 2 * sum(ifelse(y > 0, y * log(y / yhat), 0) - (y - yhat))
  crit <- data.frame(
    Criterion          = c("RMSE", "Poisson deviance (total)", "Poisson deviance / n", "MAE"),
    Poisson_GLM        = c(rmse(y24, pred_glm),  poisson_deviance(y24, pred_glm),
                           poisson_deviance(y24, pred_glm) / G,  mean(abs(y24 - pred_glm))),
    Poisson_gamma_GLMM = c(rmse(y24, pred_glmm), poisson_deviance(y24, pred_glmm),
                           poisson_deviance(y24, pred_glmm) / G, mean(abs(y24 - pred_glmm))))
  pred_tab <- data.frame(Region = d24$Region, Actual = y24,
                         GLM_pred = pred_glm, GLM_err = y24 - pred_glm,
                         GLMM_mu = mu24_glmm, E_u_given_y = u_post,
                         GLMM_pred = pred_glmm, GLMM_err = y24 - pred_glmm)
  cat(sprintf("\n=== %d hold-out predictions by region ===\n", test_year))
  print(pred_tab, digits = 3, row.names = FALSE)
  cat(sprintf("\n=== %d prediction criteria (n = %d regions) ===\n", test_year, G))
  print(crit, digits = 4, row.names = FALSE)
  cat(sprintf("Total actual = %d | GLM total = %.1f | GLMM total = %.1f\n", sum(y24), sum(pred_glm), sum(pred_glmm)))

  ## write outputs ------------------------------------------------------------
  write.csv(tab,      paste0("regression_table_", tag, ".csv"),     row.names = FALSE)
  write.csv(fit,      paste0("goodness_of_fit_", tag, ".csv"),      row.names = FALSE)
  write.csv(eb,       paste0("region_random_effects_", tag, ".csv"), row.names = FALSE)
  write.csv(pred_tab, paste0("predictions_", test_year, "_", tag, ".csv"),         row.names = FALSE)
  write.csv(crit,     paste0("prediction_criteria_", test_year, "_", tag, ".csv"), row.names = FALSE)
  invisible(list(table = tab, fit = fit, eb = eb, pred = pred_tab, crit = crit))
}

## ---- Run both specifications -----------------------------------------------
full    <- run_spec(c("RecvBudget_m1Y", "GRDPperbil", "Population", "Area"), tag = "full")
nobudget<- run_spec(c("GRDPperbil", "Population", "Area"),                   tag = "nobudget")
