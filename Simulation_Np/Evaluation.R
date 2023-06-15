# jonashaslbeck@protonmail.com; June 15th, 2023

# --------------------------------------------------------
# ------- What is happening here? ------------------------
# --------------------------------------------------------

# Analyze output of Simulation Study 2 on Permutations

# --------------------------------------------------------
# ------- Load Packages ----------------------------------
# --------------------------------------------------------

library(RColorBrewer)

# --------------------------------------------------------
# ------- Loading Sim Results ----------------------------
# --------------------------------------------------------

l_all <- list()
v_files <- list.files("Simulation_Np/output/")
n_files <- length(v_files)
l_files <- list()
for(i in 1:n_files) l_files[[i]] <- readRDS(paste0("Simulation_Np/output/", v_files[i]))

# --------------------------------------------------------
# ------- Make 2x2 Figure for Appendix -------------------
# --------------------------------------------------------

pdf("Figures/Np_simulation.pdf", width = 8, height = 7.5)

par(mfrow=c(2,2), mar=c(4.5,4.5,2.5,1))

# Pick critical values
n_crit <- 4
v_crit <- seq(0.021, 0.045, length=n_crit)
cols <- brewer.pal(n_crit, "Set1")
N <- 10000

line_tit <- 1

# ----- Panel A: Sampling Distribution -----
# Just take first iteration
X <- l_files[[1]]$SampDist$Phi_mean[2, 2, 1:N] # Note that this list entry of the output object has a different name in the latest verison of the package
plot.new()
plot.window(xlim=c(-0.07, 0.07), ylim=c(0, 30))
axis(1, seq(-0.07, 0.07, length=5))
axis(2, las=2)
title(expression(paste("Sampling Distribution with ", N[p], " = 10,000")), font.main=1, line=line_tit)
title(xlab=expression(paste("Group difference ", Delta)))
title(ylab="Density")
hist(X,
     breaks=50, xlim=c(-0.05, .05), freq = FALSE, ylim=c(0, 35), add=TRUE)
x <- seq(-.1, .1, length=N)
x_dnorm <- dnorm(x, mean=mean(X), sd=sd(X))
lines(x, x_dnorm)
abline(v=v_crit, col=cols)

# Legend
crit_round <- round(v_crit, 3)
v_leg <- LETTERS[1:4]
for(j in 1:4) {
  b_j <- bquote(Delta[crit] ~ "=" ~ .(crit_round[j]))
  text(-0.07, 30-j*2, b_j, col=cols[j], adj=0)
}

# ----- Panel B: Convergence -----
plot.new()
plot.window(xlim=c(1, N), ylim=c(0,0.2))
axis(1)
axis(2, las=2)
title("Convergence of p-values", font.main=1, line=line_tit)
title(xlab=expression(paste("Permutations ", N[p])))
title(ylab="p-value")
v_ptrue <- rep(NA, n_crit)
for(i in 1:n_crit) {
  X_ci <- cumsum(abs(X) > v_crit[i]) / (1:N)
  lines(X_ci, col=cols[i])
  v_ptrue[i] <- X_ci[N]
  abline(h=v_ptrue[i], col=cols[i], lty=2)
}
v_ptrue


# ----- Panel C: Theoretical Bounds -----
plot.new()
plot.window(xlim=c(.15,0), ylim=c(1, 15000))
axis(1, seq(.15,0, length=4))
axis(2, las=2, seq(0, 15000, length=11))
cols_con <- brewer.pal(8, "BrBG")
title("Theoretical bounds on SD(pval)", font.main=1, line=line_tit)
title(xlab="True p-value")
title(ylab=expression(paste("Permutations ", N[p])), line=3.4)

# Plot data
v_alpha <- seq(.15, 0, length=100)
v_pre <- seq(0.0025, 0.01, length=8)
for(j in 1:8) {
  v_minK <- (v_alpha*(1-v_alpha)) / v_pre[j]^2
  lines(v_alpha, v_minK, type="l", col=cols_con[j])
}

abline(v=v_ptrue, col=cols)

# Legend
# rect()
rect(.146, 6420, 0.08, 15100, col="white")
legend("topleft",
       legend=paste0("SD(pval) = ", round(v_pre, 4)),
       text.col = cols_con, bty="n")


# ----- Panel D: Theoretical vs. Empirical -----
plot.new()
plot.window(xlim=c(1, N), ylim=c(0,0.03))
axis(1)
axis(2, las=2)
title("Theoretical vs. Empirical SD(pval)", font.main=1, line=line_tit)
title(xlab=expression(paste("Permutations ", N[p])))
title(ylab="sd(pval)", line=3.4)

for(i in 1:n_crit) {
  m_X <- matrix(NA, nrow = N, ncol=10)
  for(p in 1:10) {
    X <- l_files[[p]]$SampDist$Phi_mean[1, 1, 1:N]
    X_ci <- cumsum(abs(X) > v_crit[i]) / (1:N)
    m_X[, p] <- X_ci
  }
  sd_p <- apply(m_X, 1, function(x) sd(x))
  est_true_p <- mean(m_X[N, ])
  theory_err <- sqrt((est_true_p*(1-est_true_p)) / (1:N))
  lines(sd_p, col=cols[i])
  lines(theory_err, col=cols[i], lty=2)
}

legend("topright", legend=c("Theory", "Empirical"), lty=2:1, bty="n")

dev.off()





