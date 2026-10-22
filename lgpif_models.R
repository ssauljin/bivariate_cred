## =============================================================================
## Wisconsin Local Government Property Insurance Fund (LGPIF):
## total claims N_it (all six coverages) and vehicle claims M_it (PN+PO+CN+CO).
## Estimation on 2006-2010 (data.RData), out-of-sample prediction of 2011 (dataout.RData).
##
## Total claims  : Poisson GLM vs. Poisson-gamma GLMM (gamma random effect on entity),
##                 exact marginal ML from the closed-form entity likelihood.
## Vehicle claims: binomial-probit GLM vs. binomial-probit GLMM (Gaussian random effect
##                 on entity), exact marginal ML by Gauss-Hermite quadrature.
##                 Entity-years without vehicle coverage have M_it = 0 by construction;
##                 they are excluded from the vehicle part and predicted as M = 0.
##
## Prediction for 2011:
##   N : GLM mu_i ; GLMM mu_i * E[u_i | history]
##   M : (a) given N_2011: N * p_hat_i ; (b) unconditional: E[N_2011 | history] * p_hat_i
##
## Data files are not redistributed here; place data.RData and dataout.RData
## (E. W. Frees, LGPIF) in data/lgpif/.   Base R only.   Run: Rscript lgpif_models.R
## =============================================================================

data_dir <- "data/lgpif"
out_dir  <- "results_lgpif"
dir.create(out_dir, showWarnings = FALSE)
load(file.path(data_dir, "data.RData"))      # object `data`   : 2006-2010
load(file.path(data_dir, "dataout.RData"))   # object `dataout`: 2011

prep <- function(d) {
  d$N      <- with(d, FreqBC + FreqIM + FreqPN + FreqPO + FreqCN + FreqCO)
  d$M      <- with(d, FreqPN + FreqPO + FreqCN + FreqCO)
  d$CovVeh <- with(d, CoveragePN + CoveragePO + CoverageCN + CoverageCO)
  d$HasVeh <- d$CovVeh > 0
  d
}
train <- prep(data); test <- prep(dataout)
train <- train[order(train$PolicyNum, train$Year), ]
test  <- test[order(test$PolicyNum), ]

## covariates (coverage amounts in millions); City is the baseline entity type.
## NoClaimCredit is omitted because it is itself an experience-rating variable.
design <- function(d) {
  X <- cbind(1, d$TypeCounty, d$TypeMisc, d$TypeSchool, d$TypeTown, d$TypeVillage, d$IsRC,
             log1p(d$CoverageBC), d$lnDeductBC, log1p(d$CoverageIM), log1p(d$CovVeh))
  colnames(X) <- c("(Intercept)", "TypeCounty", "TypeMisc", "TypeSchool", "TypeTown",
                   "TypeVillage", "IsRC", "log(1+CoverageBC)", "lnDeductBC",
                   "log(1+CoverageIM)", "log(1+CoverageVeh)")
  X
}

## ---- Helpers (as in binomial_probit_models.R) ------------------------------------
gh_normal <- function(Q) {
  J <- matrix(0, Q, Q)
  J[cbind(1:(Q - 1), 2:Q)] <- J[cbind(2:Q, 1:(Q - 1))] <- sqrt(1:(Q - 1))
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}
num_hessian <- function(f_grad, theta, h = 1e-5) {
  k <- length(theta); H <- matrix(0, k, k)
  for (j in seq_len(k)) { e <- rep(0, k); e[j] <- h
    H[, j] <- (f_grad(theta + e) - f_grad(theta - e)) / (2 * h) }
  (H + t(H)) / 2
}
newton_polish <- function(theta, grad_f) {
  for (it in 1:30) {
    g <- grad_f(theta); if (max(abs(g)) < 1e-8) break
    theta <- theta - solve(num_hessian(grad_f, theta), g)
  }
  theta
}
## sum the columns of a (rows x nobs) matrix within groups, returning rows x G (zeros for empty groups)
group_sum_cols <- function(A, grp, G) {
  S <- matrix(0, nrow(A), G)
  R <- rowsum(t(A), grp, reorder = TRUE)
  S[, as.integer(rownames(R))] <- t(R)
  S
}
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
mae  <- function(y, yhat) mean(abs(y - yhat))
xlogx_ratio <- function(a, b) ifelse(a > 0, a * log(a / b), 0)
pois_dev  <- function(y, yhat) 2 * sum(xlogx_ratio(y, yhat) - (y - yhat))
binom_dev <- function(m, n, p) 2 * sum(xlogx_ratio(m, n * p) + xlogx_ratio(n - m, n * (1 - p)))

ent  <- sort(unique(c(train$PolicyNum, test$PolicyNum)))
G    <- length(ent)
grpT <- match(train$PolicyNum, ent)
grpO <- match(test$PolicyNum, ent)
XT <- design(train); XO <- design(test); p <- ncol(XT)
cat(sprintf("Training: %d entity-years, %d entities with training data; test: %d entities (%d new)\n",
            nrow(train), length(unique(grpT)), nrow(test), sum(!(test$PolicyNum %in% train$PolicyNum))))
cat(sprintf("Claims in training: N = %d, M (vehicle) = %d; entity-years with vehicle coverage: %d\n",
            sum(train$N), sum(train$M), sum(train$HasVeh)))
cat(sprintf("Entity-years with M > 0 but no vehicle coverage: training %d, test %d\n",
            sum(train$M > 0 & !train$HasVeh), sum(test$M > 0 & !test$HasVeh)))

## =============================================================================
## Part 1: total claims N (Poisson GLM vs Poisson-gamma GLMM)
## =============================================================================
y <- train$N; lgam_y <- lgamma(y + 1)
gsumT <- function(x) { s <- numeric(G); r <- rowsum(x, grpT); s[as.integer(rownames(r))] <- r; s }
Y_g <- gsumT(y)
nllN <- function(theta) {
  beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
  eta <- as.vector(XT %*% beta)
  M_g <- gsumT(exp(eta))
  -(sum(y * eta - lgam_y) +
    sum(phi * logphi - lgamma(phi) + lgamma(phi + Y_g) - (phi + Y_g) * log(phi + M_g)))
}
grN <- function(theta) {
  beta <- theta[1:p]; logphi <- theta[p + 1]; phi <- exp(logphi)
  mu <- as.vector(exp(XT %*% beta))
  M_g <- gsumT(mu)
  w_g <- (phi + Y_g) / (phi + M_g)
  g_beta <- crossprod(XT, y) - crossprod(XT, mu * w_g[grpT])
  g_lphi <- phi * sum(logphi + 1 - digamma(phi) + digamma(phi + Y_g) - log(phi + M_g) - w_g)
  -c(as.vector(g_beta), g_lphi)
}
glmN <- glm(y ~ XT - 1, family = poisson(link = "log"))
names(glmN$coefficients) <- colnames(XT)
optN <- optim(c(coef(glmN), log(1)), nllN, grN, method = "BFGS",
              control = list(maxit = 10000, reltol = 1e-14))
thN <- newton_polish(optN$par, grN)
cat(sprintf("\n[N] GLMM converged | max|grad| = %.2e\n", max(abs(grN(thN)))))
VN <- solve(num_hessian(grN, thN)); seN <- sqrt(diag(VN))
betaN <- thN[1:p]; phiN <- exp(thN[p + 1])
llN_glm <- as.numeric(logLik(glmN)); llN_glmm <- -nllN(thN)
lrN <- 2 * (llN_glmm - llN_glm)
M_g <- gsumT(as.vector(exp(XT %*% betaN)))
u_post <- (phiN + Y_g) / (phiN + M_g)

## =============================================================================
## Part 2: vehicle claims M given N (probit GLM vs probit GLMM)
## =============================================================================
keep <- train$HasVeh & train$N > 0
n <- train$N[keep]; m <- train$M[keep]; X <- XT[keep, , drop = FALSE]; grp <- grpT[keep]
nobs <- sum(keep); stopifnot(all(m <= n))
cat(sprintf("\n[M] binomial part: %d entity-years with vehicle coverage and N > 0, %d entities, %d claims, %d vehicle\n",
            nobs, length(unique(grp)), sum(n), sum(m)))
glmM <- glm(cbind(m, n - m) ~ X - 1, family = binomial(link = "probit"))
names(glmM$coefficients) <- colnames(X)
pearsonM <- sum(residuals(glmM, type = "pearson")^2) / glmM$df.residual

Q_fit <- 200
gh <- gh_normal(Q_fit); lchoose_nm <- sum(lchoose(n, m))
node_terms <- function(beta, sigma) {
  eta <- as.vector(X %*% beta)
  A <- outer(sigma * gh$z, eta, "+")
  list(A = A, ll = sweep(pnorm(A, log.p = TRUE), 2, m, "*") +
                   sweep(pnorm(-A, log.p = TRUE), 2, n - m, "*"))
}
region_logpost <- function(ll) sweep(group_sum_cols(ll, grp, G), 1, log(gh$w), "+")
logsumexp_cols <- function(L) { mx <- apply(L, 2, max); mx + log(colSums(exp(sweep(L, 2, mx)))) }
used <- sort(unique(grp))                               # entities contributing to the likelihood
nllM <- function(theta) {
  nt <- node_terms(theta[1:p], exp(theta[p + 1]))
  -(lchoose_nm + sum(logsumexp_cols(region_logpost(nt$ll)[, used, drop = FALSE])))
}
grM <- function(theta) {
  beta <- theta[1:p]; sigma <- exp(theta[p + 1])
  nt <- node_terms(beta, sigma); L <- region_logpost(nt$ll)
  P <- exp(sweep(L, 2, logsumexp_cols(L)))
  mills <- function(a) exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
  D <- sweep(mills(nt$A), 2, m, "*") - sweep(mills(-nt$A), 2, n - m, "*")
  Pobs <- P[, grp, drop = FALSE]
  -c(as.vector(crossprod(X, colSums(Pobs * D))), sigma * sum(colSums(Pobs * D * gh$z)))
}
optM <- optim(c(coef(glmM), log(0.5)), nllM, grM, method = "BFGS",
              control = list(maxit = 10000, reltol = 1e-14))
thM <- newton_polish(optM$par, grM)
cat(sprintf("[M] GLMM converged | max|grad| = %.2e\n", max(abs(grM(thM)))))
VM <- solve(num_hessian(grM, thM)); seM <- sqrt(diag(VM))
betaM <- thM[1:p]; sigM <- exp(thM[p + 1])
llM_glm <- as.numeric(logLik(glmM)); llM_glmm <- -nllM(thM)
lrM <- 2 * (llM_glmm - llM_glm)
ll_check <- { gh_save <- gh; gh <- gh_normal(400); v <- -nllM(thM); gh <- gh_save; v }
eta_hat <- as.vector(X %*% betaM)
ll_int <- lchoose_nm + sum(sapply(used, function(g) {
  idx <- grp == g
  lf <- function(z) sum(m[idx] * pnorm(eta_hat[idx] + sigM * z, log.p = TRUE) +
                        (n[idx] - m[idx]) * pnorm(-eta_hat[idx] - sigM * z, log.p = TRUE)) + dnorm(z, log = TRUE)
  zs <- seq(-8, 8, by = 0.01); lv <- sapply(zs, lf); mx <- max(lv); z0 <- zs[which.max(lv)]
  f <- function(z) exp(sapply(z, lf) - mx)
  mx + log(integrate(f, max(-10, z0 - 3), min(10, z0 + 3), rel.tol = 1e-10, subdivisions = 2000L)$value +
           integrate(f, -10, max(-10, z0 - 3), rel.tol = 1e-10, subdivisions = 2000L)$value +
           integrate(f, min(10, z0 + 3), 10, rel.tol = 1e-10, subdivisions = 2000L)$value)
}))
cat(sprintf("[M] loglik GH(%d) = %.6f | GH(400) = %.6f | integrate() = %.6f\n", Q_fit, llM_glmm, ll_check, ll_int))

## posterior predictive probability of a vehicle claim (entities without history keep the prior)
gh <- gh_normal(400)
L <- region_logpost(node_terms(betaM, sigM)$ll)
P <- exp(sweep(L, 2, logsumexp_cols(L)))                 # Q x G; prior weights for empty groups
vq <- sigM * gh$z
v_post <- as.vector(crossprod(P, vq))

## =============================================================================
## Tables
## =============================================================================
coef_table <- function(glm_fit, theta, se, re_name, re_est, re_se) {
  k <- length(coef(glm_fit)); gse <- sqrt(diag(vcov(glm_fit)))
  tab <- data.frame(Term = names(coef(glm_fit)),
                    GLM_Estimate = coef(glm_fit), GLM_SE = gse, GLM_p = 2 * pnorm(-abs(coef(glm_fit) / gse)),
                    GLMM_Estimate = theta[1:k], GLMM_SE = se[1:k],
                    GLMM_p = 2 * pnorm(-abs(theta[1:k] / se[1:k])), row.names = NULL)
  rbind(tab, data.frame(Term = re_name, GLM_Estimate = NA, GLM_SE = NA, GLM_p = NA,
                        GLMM_Estimate = re_est, GLMM_SE = re_se, GLMM_p = NA))
}
tabN <- coef_table(glmN, thN, seN, "Var(u) = 1/phi", 1 / phiN, seN[p + 1] / phiN)
tabM <- coef_table(glmM, thM, seM, "Var(v) = sigma^2", sigM^2, 2 * sigM^2 * seM[p + 1])
nN <- nrow(train)
fit <- data.frame(
  Statistic = c("Observations", "Parameters", "Log-likelihood", "AIC", "BIC", "LR test", "LR p (0.5 chi2_1)",
                "Pearson chi2/df (GLM)"),
  N_GLM  = c(nN, p, llN_glm, -2 * llN_glm + 2 * p, -2 * llN_glm + log(nN) * p, NA, NA,
             sum(residuals(glmN, type = "pearson")^2) / glmN$df.residual),
  N_GLMM = c(nN, p + 1, llN_glmm, -2 * llN_glmm + 2 * (p + 1), -2 * llN_glmm + log(nN) * (p + 1),
             lrN, 0.5 * pchisq(lrN, 1, lower.tail = FALSE), NA),
  M_GLM  = c(nobs, p, llM_glm, -2 * llM_glm + 2 * p, -2 * llM_glm + log(nobs) * p, NA, NA, pearsonM),
  M_GLMM = c(nobs, p + 1, llM_glmm, -2 * llM_glmm + 2 * (p + 1), -2 * llM_glmm + log(nobs) * (p + 1),
             lrM, 0.5 * pchisq(lrM, 1, lower.tail = FALSE), NA))
cat("\n=== Total claims N: coefficients ===\n");   print(tabN, digits = 4, row.names = FALSE)
cat("\n=== Vehicle claims M: coefficients ===\n"); print(tabM, digits = 4, row.names = FALSE)
cat("\n=== Fit statistics ===\n");                 print(fit, digits = 6, row.names = FALSE)

## =============================================================================
## 2011 predictions
## =============================================================================
N11 <- test$N; M11 <- test$M; veh <- test$HasVeh
EN_glm  <- as.vector(exp(XO %*% coef(glmN)))
EN_glmm <- as.vector(exp(XO %*% betaN)) * u_post[grpO]
eta11   <- as.vector(XO %*% betaM)
p_glm   <- ifelse(veh, as.vector(pnorm(XO %*% coef(glmM))), 0)
p_glmm  <- ifelse(veh, colSums(P[, grpO, drop = FALSE] * pnorm(outer(vq, eta11, "+"))), 0)

critN <- data.frame(Criterion = c("RMSE", "MAE", "Poisson deviance", "Predicted total"),
                    GLM  = c(rmse(N11, EN_glm), mae(N11, EN_glm), pois_dev(N11, EN_glm), sum(EN_glm)),
                    GLMM = c(rmse(N11, EN_glmm), mae(N11, EN_glmm), pois_dev(N11, EN_glmm), sum(EN_glmm)))
ok <- veh & N11 > 0
critM <- data.frame(
  Criterion = c("(a) given N: RMSE", "(a) given N: MAE", "(a) given N: binomial deviance*",
                "(b) unconditional: RMSE", "(b) unconditional: MAE", "(b) unconditional: Poisson deviance",
                "Predicted total, given N", "Predicted total, unconditional"),
  GLM  = c(rmse(M11, N11 * p_glm), mae(M11, N11 * p_glm), binom_dev(M11[ok], N11[ok], p_glm[ok]),
           rmse(M11, EN_glm * p_glm), mae(M11, EN_glm * p_glm), pois_dev(M11[veh], (EN_glm * p_glm)[veh]),
           sum(N11 * p_glm), sum(EN_glm * p_glm)),
  GLMM = c(rmse(M11, N11 * p_glmm), mae(M11, N11 * p_glmm), binom_dev(M11[ok], N11[ok], p_glmm[ok]),
           rmse(M11, EN_glmm * p_glmm), mae(M11, EN_glmm * p_glmm), pois_dev(M11[veh], (EN_glmm * p_glmm)[veh]),
           sum(N11 * p_glmm), sum(EN_glmm * p_glmm)))
cat(sprintf("\n=== 2011 prediction: total claims N (%d entities, actual total %d) ===\n", nrow(test), sum(N11)))
print(critN, digits = 5, row.names = FALSE)
cat(sprintf("\n=== 2011 prediction: vehicle claims M (actual total %d; %d entities with vehicle coverage) ===\n",
            sum(M11), sum(veh)))
cat("* binomial deviance over entities with vehicle coverage and N > 0; Poisson deviance over entities with vehicle coverage\n")
print(critM, digits = 5, row.names = FALSE)

## credibility spread
cat(sprintf("\nE[u|H] over test entities: %s\n", paste(sprintf("%.3f", quantile(u_post[grpO], c(0, .1, .5, .9, 1))), collapse = " / ")))
cat(sprintf("E[v|H] over test entities with vehicle coverage: %s\n",
            paste(sprintf("%.3f", quantile(v_post[grpO][veh], c(0, .1, .5, .9, 1))), collapse = " / ")))

pred <- data.frame(PolicyNum = test$PolicyNum, HasVeh = veh, N_2011 = N11, M_2011 = M11,
                   EN_GLM = EN_glm, EN_GLMM = EN_glmm, E_u_given_H = u_post[grpO], E_v_given_H = v_post[grpO],
                   p_GLM = p_glm, p_GLMM = p_glmm,
                   M_given_N_GLM = N11 * p_glm, M_given_N_GLMM = N11 * p_glmm,
                   M_uncond_GLM = EN_glm * p_glm, M_uncond_GLMM = EN_glmm * p_glmm)
write.csv(tabN, file.path(out_dir, "regression_table_N.csv"), row.names = FALSE)
write.csv(tabM, file.path(out_dir, "regression_table_M.csv"), row.names = FALSE)
write.csv(fit,  file.path(out_dir, "goodness_of_fit.csv"),    row.names = FALSE)
write.csv(critN, file.path(out_dir, "prediction_criteria_N_2011.csv"), row.names = FALSE)
write.csv(critM, file.path(out_dir, "prediction_criteria_M_2011.csv"), row.names = FALSE)
write.csv(pred, file.path(out_dir, "predictions_2011.csv"), row.names = FALSE)

## ---- Robustness of the 2011 comparison --------------------------------------------
big <- N11 >= quantile(N11, 0.99); s <- !big; v <- veh & s
cat(sprintf("\n=== Robustness: excluding the %d entities in the top 1%% of N_2011 (N >= %d) ===\n",
            sum(big), min(N11[big])))
cat(sprintf("N RMSE: GLM %.3f | GLMM %.3f\n", rmse(N11[s], EN_glm[s]), rmse(N11[s], EN_glmm[s])))
cat(sprintf("M RMSE (a): GLM %.3f | GLMM %.3f;  (b): GLM %.3f | GLMM %.3f\n",
            rmse(M11[v], (N11 * p_glm)[v]), rmse(M11[v], (N11 * p_glmm)[v]),
            rmse(M11[v], (EN_glm * p_glm)[v]), rmse(M11[v], (EN_glmm * p_glmm)[v])))
set.seed(1)
boot <- replicate(2000, {
  i <- sample(nrow(test), replace = TRUE); vv <- veh[i]
  c(rmse(N11[i], EN_glm[i]) - rmse(N11[i], EN_glmm[i]),
    rmse(M11[i][vv], (N11 * p_glm)[i][vv]) - rmse(M11[i][vv], (N11 * p_glmm)[i][vv]),
    rmse(M11[i][vv], (EN_glm * p_glm)[i][vv]) - rmse(M11[i][vv], (EN_glmm * p_glmm)[i][vv]))
})
cat(sprintf("Entity bootstrap (2000 resamples), P(RMSE_GLMM < RMSE_GLM): N %.3f | M (a) %.3f | M (b) %.3f\n",
            mean(boot[1, ] > 0), mean(boot[2, ] > 0), mean(boot[3, ] > 0)))

## ---- Descriptive: claims by year ------------------------------------------------------
allyrs <- rbind(train[, c("Year", "N", "M")], test[, c("Year", "N", "M")])
by_year <- aggregate(cbind(N, M) ~ Year, data = allyrs, FUN = sum)
by_year$Entities <- as.vector(table(allyrs$Year)); by_year$Share <- by_year$M / by_year$N
cat("\n=== Claims by year ===\n"); print(by_year, digits = 3, row.names = FALSE)
write.csv(by_year, file.path(out_dir, "claims_by_year.csv"), row.names = FALSE)
