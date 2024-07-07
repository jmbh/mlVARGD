# jonashaslbeck@protonmail.com; March 2nd, 2023

# -------------------------------------------------
# -------- What is happening here? ----------------
# -------------------------------------------------

# Some auxiliary functions



# -------------------------------------------------
# -------- Estimate VAR ---------------------------
# -------------------------------------------------

EstimateVAR <- function(data, roundpoints = 4, summaries = FALSE, changescore = FALSE) {

  p <- ncol(data)
  Phi <- matrix(NA, p, p)
  alpha <- rep(NA, p)
  residuals <- matrix(NA, nrow(data)-1, p)
  summary_list<- list()

  for(i in 1:p) {

    # Estimates the change-score version of the var model instead
    if(changescore == TRUE){
      y = data[-1,i] - data[-nrow(data), i]
    } else {
      y <- data[-1, i]
    }


    # predicted cases
    X <- data[-nrow(data), ]

    coefs <- lm(y ~ X)
    Phi[i, ] <- coefs$coefficients[-1]
    alpha[i] <- coefs$coefficients[1]
    residuals[,i] <- coefs$residuals
    if(summaries) summary_list[[i]] <- summary(coefs)
  }

  Psi <- cov(residuals)
  mu <- as.vector(solve(diag(p)-Phi)%*%alpha)

  coefs <- list(round(alpha, roundpoints),
                round(mu,roundpoints),
                round(Phi, roundpoints),
                round(Psi,roundpoints))

  names(coefs) <- c("intercepts", "means","phi","psi")

  if(summaries){
    return(summary_list)
  } else {
    return(coefs)
  }

} # eoF

