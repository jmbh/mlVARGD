# jonashaslbeck@protonmail.com; July 7th, 2024

# --------------------------------------------------------
# ------- What is happening here? ------------------------
# --------------------------------------------------------

# Analyze output of Simulation Study 1 on Permutations

# --------------------------------------------------------
# ------- Load Packages ----------------------------------
# --------------------------------------------------------

library(scales)
library(RColorBrewer)
library(reshape2)
library(xtable)

# --------------------------------------------------------
# ------- Some Aux for Evaluation ------------------------
# --------------------------------------------------------

plotLabel <- function(x, srt=0, col="black",
                      xpos=.6, ypos=.5, cex=1.4) {
  par(mar=rep(0, 4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(xpos, ypos, x, srt=srt, cex=cex, col=col)
}

# --------------------------------------------------------
# ------- Loading Sim Results ----------------------------
# --------------------------------------------------------

l_all <- list()
v_files <- list.files(paste0("Simulation_Main_New/output/"))
n_files <- length(v_files)
l_files <- list()
for(i in 1:n_files) l_files[[i]] <- readRDS(paste0("Simulation_Main_New/output/", v_files[i]))

date_of_analysis <- "July1st_24"


# --------------------------------------------------------
# ------- Some Sanity Checks -----------------------------
# --------------------------------------------------------


# --------------------------------------------------------
# ------- Global Settings --------------------------------
# --------------------------------------------------------

# --- Sequence of evaluated p-values ---
v_pvals <- seq(0, 1, length=100)
# Add the three standard cutoffs
v_pvals <- c(v_pvals[1], 0.01, v_pvals[2:5], 0.05, v_pvals[6:10], 0.1, v_pvals[11:100])
n_pvals <- length(v_pvals)
v_pval_loc <- c(which(v_pvals==0.01), which(v_pvals==0.05), which(v_pvals==0.10))

v_delta <- c(0.04, 0.08, 0.12)

# --------------------------------------------------------
# ------- Preprocess Results: Parametric Test ------------
# --------------------------------------------------------

nIter <- n_files
n_delta <- 3
n_Nsubj <- 3
n_prop <- 2


## Storage
# iterations, n_subj, n_prop, n_delta, 3 para types, 100 x pval thresh, 2 outcomes (TPR, FPR)
a_res_para <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta, 3, n_pvals, 2))
a_bias_VAR_diff <- array(NA, dim=c(8, 8, nIter, n_Nsubj, n_prop, n_delta, 2)) # Last dim: True vs. Estimated
a_resamples1 <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta)) # Resamples in data generation (to achieve sample VAR models)
a_resamples2 <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta)) # Resamples in data generation (to achieve pos def residual covariances)

for(i in 1:nIter) {
  j <- 1   # counter for design matrix
  for(i_d in 1:3) {
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {
        for(type in c(1,2,3)) {
          for(i_pval in 1:n_pvals) {

            # Processing conditional on structure of matric/no of parameters
            if(type %in% 1) {
              diff_true <- l_files[[i]]$datagen[[j]]$true$GroupDiffs[[type]] != 0
              diff_est <- l_files[[i]]$outcomes_para[[j]]$Pval[[type]] < v_pvals[i_pval]
            } else {
              diff_true <- l_files[[i]]$datagen[[j]]$true$GroupDiffs[[c(1,3,5)[type]]]
              diff_true <- diff_true[lower.tri(diff_true)] != 0
              diff_est <- l_files[[i]]$outcomes_para[[j]]$Pval[[type]]
              diff_est <- diff_est[lower.tri(diff_est)] < v_pvals[i_pval]
            }

            # Get measures
            TPR <- mean(diff_est[diff_true])
            FPR <- mean(diff_est[!diff_true])


            # Save
            a_res_para[i, i_sub, i_pr, i_d, type, i_pval, 1:2] <- c(TPR, FPR)

          }
        }

        # Save no. of resamples
        a_resamples1[i, i_sub, i_pr, i_d] <- l_files[[i]]$datagen[[j]]$dataSamp$counter_base
        a_resamples2[i, i_sub, i_pr, i_d] <- l_files[[i]]$datagen[[j]]$dataSamp$counter_datasamp

        ## Save estimate of group differences in lagged effects (to investigate alpha-misalignment)
        # True
        a_bias_VAR_diff[, , i, i_sub, i_pr, i_d, 1] <- l_files[[i]]$datagen[[j]]$true$GroupDiffs$m_GD_phi
        # Est
        a_bias_VAR_diff[, , i, i_sub, i_pr, i_d, 2] <- l_files[[i]]$outcomes_para[[j]]$EmpDiffs$Lagged_fixed

        j <- j+1 # update counter; should be 18+1 in the end

        print(paste0("i= ", i, "; j=", j, "; i_sub=", i_sub, "; i_pr=", i_pr, "; i_d=", i_d))

      }
    }
  }
} # end for: iters

# Aggregate across repetitions of simulation study
a_res_para_agg <- apply(a_res_para, 2:7, mean)
dim(a_res_para)



# --------------------------------------------------------
# ------- Preprocess Results: Permutation Test -----------
# --------------------------------------------------------

nIter <- n_files
n_delta <- 3
n_Nsubj <- 3
n_prop <- 2

## Storage
# iterations, n_subj, n_prop, n_delta, 5 vartypes, 100 x pval thresh, 2 outcomes (TPR, FPR)
a_res <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta, 5, n_pvals, 2)) # Results

for(i in 1:nIter) {
  j <- 1   # counter for design matrix
  for(i_d in 1:3) { # delta
    for(i_pr in 1:2) { # proportion
      for(i_sub in 1:3) { # number of subjects


        for(type in 1:5) {
          for(i_pval in 1:n_pvals) {

            # Processing conditional on structure of matric/no of parameters
            if(type %in% 1:2) {
              diff_true <- l_files[[i]]$datagen[[j]]$true$GroupDiffs[[type]] != 0
              diff_est <- l_files[[i]]$outcomes[[j]]$Pval[[type]] < v_pvals[i_pval]
            } else {
              diff_true <- l_files[[i]]$datagen[[j]]$true$GroupDiffs[[type]]
              diff_true <- diff_true[lower.tri(diff_true)] != 0
              diff_est <- l_files[[i]]$outcomes[[j]]$Pval[[type]]
              diff_est <- diff_est[lower.tri(diff_est)] < v_pvals[i_pval]
            }

            # Get measures
            TPR <- mean(diff_est[diff_true])
            FPR <- mean(diff_est[!diff_true])

            # Save
            a_res[i, i_sub, i_pr, i_d, type, i_pval, 1:2] <- c(TPR, FPR)

          }
        }

        print(paste0("i= ", i, "; j=", j, "; i_sub=", i_sub, "; i_pr=", i_pr, "; i_d=", i_d))

        j <- j+1 # update counter; should be 18+1 in the end

      }
    }
  }
} # end for: iter

# Aggregate across repetitions of simulation study
a_res_agg <- apply(a_res, 2:7, mean)



# --------------------------------------------------------
# ------- Main Results Figure: ROC Curves Parametric -----
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

pdf(paste0("Figures/Fig_Res_Parametric_", date_of_analysis,".pdf"), width=7, height=9*(3/5))

# ----- Make Layout -----
lmat <- matrix(7:15, 3, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:6), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 3)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Fixed: Lagged", srt=90, ypos=.6)
plotLabel("Fixed: Innovation", srt=90, ypos=.6)
plotLabel("Between", srt=90, ypos=.6)

# ----- Plot Data -----

for(type in 1:3) { # type of parameter
  for(i_d in 1:3) { # delta

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="FPR", line=2.1)
    title(ylab="TPR", line=2.3)
    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {

        lines(a_res_para_agg[i_sub, i_pr, i_d, type, , 2], a_res_para_agg[i_sub, i_pr, i_d, type, , 1],
              col=scales::alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)

        # Add alphas for 0.01, 0.05, 0.10
        for(alpha in 1:3) points(a_res_para_agg[i_sub, i_pr, i_d, type, v_pval_loc[alpha], 2],
                                 a_res_para_agg[i_sub, i_pr, i_d, type, v_pval_loc[alpha], 1],
                                 col = cols[i_sub],
                                 pch=alpha+3)


      }
    }
    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.6, .29, legend=c("50:50", "20:80"), bty="n", lty=1:2, cex=cex_size)
      legend(.25, .35, legend=c(expression(paste(alpha==0.01)),
                                expression(paste(alpha==0.05)),
                                expression(paste(alpha==0.10))), bty="n", pch=4:6, cex=cex_size)
    }
  }
}
dev.off()

# ----- Get some specific values to report in paper -----
alpha <- 2 # alpha = 0.05
i_d <- 1 # Delta
type <- 1 # type of parameter
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2) # TPR
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2) # FPR
i_d <- 2
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2) # TPR
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2) # FPR
i_d <- 3
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2) # TPR
round(a_res_para_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2) # FPR



# --------------------------------------------------------
# ------- Main Results Figure: ROC Curves Permutation ----
# --------------------------------------------------------

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")


pdf(paste0("Figures/Fig_Res_Perm_", date_of_analysis, ".pdf"), width=7, height=9)

# ----- Make Layout -----
lmat <- matrix(9:23, 5, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:8), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 4)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Fixed: Lagged", srt=90, ypos=.6)
plotLabel("RE: Lagged", srt=90, ypos=.6)
plotLabel("Fixed: Innovation", srt=90, ypos=.6)
plotLabel("RE: Innovation", srt=90, ypos=.6)
plotLabel("Between", srt=90, ypos=.6)

# ----- Plot Data -----

for(type in 1:5) {
  for(i_d in 1:3) {

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="FPR", line=2.1)
    title(ylab="TPR", line=2.3)
    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {
        lines(a_res_agg[i_sub, i_pr, i_d, type, , 2], a_res_agg[i_sub, i_pr, i_d, type, , 1],
              col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)

        # Add alphas for 0.01, 0.05, 0.10
        for(alpha in 1:3) points(a_res_agg[i_sub, i_pr, i_d, type, v_pval_loc[alpha], 2],
                                 a_res_agg[i_sub, i_pr, i_d, type, v_pval_loc[alpha], 1],
                                 col = cols[i_sub],
                                 pch = alpha)
      }
    }
    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.6, .29, legend=c("50:50", "20:80"), bty="n", lty=1:2, cex=cex_size)
      legend(.25, .35, legend=c(expression(paste(alpha==0.01)),
                                expression(paste(alpha==0.05)),
                                expression(paste(alpha==0.10))), bty="n", pch=1:3, cex=cex_size)
    }
  }
}
dev.off()

# Get some specific values to report in paper
alpha <- 2 # alpha = 0.05
i_d <- 1
i_sub <- 1
type <- 1
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2)
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2)
i_d <- 2
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2)
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2)
i_d <- 3
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 1], 2)
round(a_res_agg[1:3, i_pr, i_d, type, v_pval_loc[alpha], 2], 2)


# --------------------------------------------------------
# ------- Figure: Difference Between Para vs. Non-Para ---
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

# Get subset of permutation results
a_res_nonpara_subset <- a_res_agg[, , , c(1,3,5), , ]

pdf(paste0("Figures/Fig_Res_Para_vs_NonPara_",date_of_analysis,".pdf"), width=7, height=9*(2/5))

# ----- Make Layout -----
lmat <- matrix(6:11, 2, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:5), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 2)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Fixed: Lagged", srt=90, ypos=.6)
plotLabel("Fixed: Innovation", srt=90, ypos=.6)
# plotLabel("Between", srt=90)

# ----- Plot Data -----

for(type in 1:2) { # type of parameter
  for(i_d in 1:3) { # delta

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="FPR", line=2.1)
    title(ylab="TPR", line=2.3)
    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {

      if(i_pr==1) a_res_disp <- a_res_nonpara_subset # results permutation test
      if(i_pr==2) a_res_disp <- a_res_para_agg # results parametric test

      # Average across balanced & unbalanced
      a_res_disp_aob <- apply(a_res_disp, c(1,3:6), mean)

      # balance <- 2

      for(i_sub in 1:3) {

        lines(a_res_disp_aob[i_sub, i_d, type, , 2], a_res_disp_aob[i_sub, i_d, type, , 1],
              col=scales::alpha(cols[i_sub], alpha=0.5), lty=c(1,4)[i_pr], lwd=2)


        # Add alphas for 0.01, 0.05, 0.10
        if(i_pr==1) extra <- 0
        if(i_pr==2) extra <- 3

        for(alpha in 1:3) points(a_res_disp_aob[i_sub, i_d, type, c(2,6,11)[alpha], 2],
                                 a_res_disp_aob[i_sub, i_d, type, c(2,6,11)[alpha], 1],
                                 col = cols[i_sub],
                                 pch = alpha+extra)
      }
    }

    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.4, .3, legend=c("Non-Parametric", "Parametric"), bty="n", lty=c(1,4), cex=0.8)
    }

    if(type==1 & i_d==2) {

      rect(.13, 0, .9, .38, border=FALSE, col="white")
      cex_size <- 0.9
      text(.76, .34, "Non-parametric:", cex=cex_size, adj=0.5)
      legend(.5, .35, legend=c(expression(paste(alpha==0.01)),
                               expression(paste(alpha==0.05)),
                               expression(paste(alpha==0.10))), bty="n", pch=1:3, cex=cex_size)

      text(.31, .35, "Parametric:", cex=cex_size, adj=0.5)
      legend(.13, .35, legend=c(expression(paste(alpha==0.01)),
                                expression(paste(alpha==0.05)),
                                expression(paste(alpha==0.10))), bty="n", pch=4:6, cex=cex_size)
    }

  }
}

dev.off()


# --------------------------------------------------------
# ------ Specified vs. empirical alpha/FPR [parametric] --
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

pdf(paste0("Figures/Fig_AlphaEval_Para_", date_of_analysis,".pdf"), width=7, height=9*(2/5))

# ----- Make Layout -----
lmat <- matrix(6:(6+6-1), 2, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:5), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 2)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Fixed: Lagged", srt=90, ypos=.6)
plotLabel("Fixed: Innovation", srt=90, ypos=.6)
# plotLabel("Between", srt=90)

# ----- Plot Data -----

for(type in 1:2) { # type of parameter
  for(i_d in 1:3) { # delta

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="Specified alpha", line=2.1)
    title(ylab="Observed FPR", line=2.3)
    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {

        lines(v_pvals, a_res_para_agg[i_sub, i_pr, i_d, type, , 2],
              col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)

      }
    }
    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.6, .29, legend=c("50:50", "20:80"), bty="n", lty=1:2, cex=cex_size)
      # legend(.25, .35, legend=c(expression(paste(alpha==0.01)),
      #                           expression(paste(alpha==0.05)),
      #                           expression(paste(alpha==0.10))), bty="n", pch=1:3, cex=cex_size)
    }
  }
}
dev.off()


# i_sub <- 3
# i_pr <- 1
# i_d <- 3
# type <- 1
# a_res_para_agg[i_sub, i_pr, i_d, type, 6, 2]



# --------------------------------------------------------
# ------ Specified vs. empirical alpha/FPR [permutation] -
# --------------------------------------------------------

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")


pdf(paste0("Figures/Fig_AlphaEval_NonPara_", date_of_analysis,".pdf"), width=7, height=9*(4/5))

# ----- Make Layout -----
lmat <- matrix(8:(8+3*4-1), 4, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:7), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 3)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Fixed: Lagged", srt=90, ypos=.6)
plotLabel("RE: Lagged", srt=90, ypos=.6)
plotLabel("Fixed: Innovation", srt=90, ypos=.6)
plotLabel("RE: Innovation", srt=90, ypos=.6)
# plotLabel("Between", srt=90)

# ----- Plot Data -----

for(type in 1:4) { # no between network
  for(i_d in 1:3) {

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="Specified alpha", line=2.1)
    title(ylab="Observed FPR", line=2.3)
    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {
        # lines(v_pvals, a_res_agg[i_sub, i_pr, i_d, type, , 1],
        #       col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)
        lines(v_pvals, a_res_agg[i_sub, i_pr, i_d, type, , 2], # FPR
              col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)
      }
    }
    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.6, .29, legend=c("50:50", "20:80"), bty="n", lty=1:2, cex=cex_size)
      # legend(.25, .35, legend=c(expression(paste(alpha==0.01)),
      #                           expression(paste(alpha==0.05)),
      #                           expression(paste(alpha==0.10))), bty="n", pch=1:3, cex=cex_size)
    }
  }
}
dev.off()


# --------------------------------------------------------------------
# ------- Analyze Sampl. Dists of Diffs in Fixed Lagged effects [PAPER VERSION] ------
# --------------------------------------------------------------------


cols_SDs <- brewer.pal(4, "Set2")

pdf(paste0("Figures/SIMDATA_SampDists_PAPERV_", date_of_analysis, ".pdf"), width=9, height=9)

# ----- Make Layout -----
lmat <- matrix(7:15, 3, 3, byrow = TRUE)
lmat <- rbind(1:3, lmat)
lmat <- cbind(c(0, 4:6), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 3)),
             heights = c(.09, rep(1, 3)))

# ----- Plot Labels -----

# Top row
plotLabel(expression(paste(Delta, " = 0.04")))
plotLabel(expression(paste(Delta, " = 0.08")))
plotLabel(expression(paste(Delta, " = 0.12")))
# Left Col
plotLabel("Nsubj = 50", srt=90, ypos=.6)
plotLabel("Nsubj = 100", srt=90, ypos=.6)
plotLabel("Nsubj = 200", srt=90, ypos=.6)

par(mar=c(4,3,2,1))

for(i_sub in 1:3) {
  for(i_d in 1:3) {
    phi_true <- a_bias_VAR_diff[, , , i_sub, 1, i_d, 1] # i_pr <- 1 FIXED
    phi_est <- a_bias_VAR_diff[, , , i_sub, 1, i_d, 2] # i_pr <- 1 FIXED

    # times -1, because we have diff=G1-G2 in the package but diff=G2-G1 in the sim function
    a_diff <- phi_true - (-1 * phi_est)

    # SUBSETTING
    a_diff[phi_true!=0] <- NA

    v_diff <- as.numeric(a_diff)

    y_ax_mp <- 1.3
    sgm_max <- 500

    # n_dec <- 4

    hist(v_diff, xlim=c(-.2, .2),
         breaks=seq(-.7, .7, length=101),
         xlab="Difference: Group 1 - Group 2", main="", ylim=c(0,800))

    mn <- mean(v_diff, na.rm=TRUE)
    segments(x0 = mn, y0=0, x1 = mn, y1=sgm_max, col=cols_SDs[1], lwd=3)
    text(-0.2, 700*y_ax_mp-100, paste0("Mean = ", round(mn, 6)), adj=0, col=cols_SDs[1])
    # Quantiles
    qts <- round(quantile(v_diff, probs = c(0.05, 0.95), na.rm=TRUE), 3)
    text(-0.2, 650*y_ax_mp-100, paste0("5% Quantile = ", qts[1]), adj=0, col=cols_SDs[2])
    text(-0.2, 600*y_ax_mp-100, paste0("95% Quantile = ", qts[2]), adj=0, col=cols_SDs[3])
    text(-0.2, 550*y_ax_mp-100, paste0("SD = ", round(sd(v_diff, na.rm=TRUE), 6)), col=cols_SDs[4],  adj=0)
    segments(x0 = qts[1], y0=0, x1 = qts[1], y1=sgm_max, col=cols_SDs[2], lwd=3)
    segments(x0 = qts[2], y0=0, x1 = qts[2], y1=sgm_max, col=cols_SDs[3], lwd=3)

  }
}

dev.off()


# --------------------------------------------------------
# ------- Figure to explain weird RE: Residual Results ------
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

sc <- 1.3
pdf(paste0("Figures/Fig_WeirdResult_", date_of_analysis,".pdf"), width=sc*7*(2/3)*0.9, height=sc*9*(1/5))

# ----- Make Layout -----
lmat <- matrix(4:5, 1, 2, byrow = TRUE)
lmat <- rbind(1:2, lmat)
lmat <- cbind(c(0, 3), lmat)

lo <- layout(lmat,
             widths = c(.09, rep(1, 2)),
             heights = c(.09, rep(1, 1)))
# layout.show(lo)

# ----- Plot Labels -----

# Top row
plotLabel("Alpha vs. FPR")
plotLabel("Alpha vs. TPR")
# Left Col
plotLabel("RE: Innovation", srt=90, ypos=.6)


# ----- Plot Data -----

# Display either alpha X FPR (1) or alpha X TPR (3)
for(disp in 1:2) {

  for(type in 4) { # type of parameter
    # for(i_d in 1:3) { # delta

    i_d <- 2 # fixed

    # Canvas
    par(mar=c(3.5,3.5,0.5,1))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    axis(1)
    axis(2, las=2)
    title(xlab="Specified alpha", line=2.1)

    if(disp==1) title(ylab="Observed FPR", line=2.3)
    if(disp==2) title(ylab="Observed TPR", line=2.3)

    segments(0,0,1,1, lwd=2, col="grey")

    abline(h=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)
    abline(v=seq(0, 1, length=6), lty=2, col="grey", lwd=0.5)

    # Data
    for(i_pr in 1:2) {
      for(i_sub in 1:3) {

        if(disp==1) lines(v_pvals, a_res_agg[i_sub, i_pr, i_d, type, , 2],
                          col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)
        if(disp==2) lines(v_pvals, a_res_agg[i_sub, i_pr, i_d, type, , 1],
                          col=alpha(cols[i_sub], alpha=0.5), lty=i_pr, lwd=2)

      }
    }
    # Legends
    if(disp==1 & i_d==2) {
      cex_size <- 0.9
      legend("topleft", legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend("left", legend=c("50:50", "20:80"), bty="n", lty=1:2, cex=cex_size)
    }
  }
}


dev.off()


# --------------------------------------------------------
# ------- Evaluate SE of Estimates -----------------------
# --------------------------------------------------------

# ------ Parametric Test ------
a_res_para_Big_agg_sd <- apply(a_res_para[, , , , , 7, ], c(2,5,6), function(x) sd(x, na.rm=TRUE))
dimnames(a_res_para_Big_agg_sd) <- list(c(50, 100, 200),
                                        c("Fixed: Lagged", "Fixed: Innovation", "Between"),
                                        c("TPR", "FPR"))

a_res_para_Big_agg_sd_long <- melt(a_res_para_Big_agg_sd)
colnames(a_res_para_Big_agg_sd_long) <- c("#Subjects","ParType", "Outcome", "SD")
n <- 50
a_res_para_Big_agg_sd_long$SE <- a_res_para_Big_agg_sd_long$SD / sqrt(n)
a_res_para_Big_agg_sd_long$SE <- round(a_res_para_Big_agg_sd_long$SE, 3)
a_res_para_Big_agg_sd_long$SD <- round(a_res_para_Big_agg_sd_long$SD, 3)

# Print the LaTeX table
latex_table <- xtable(a_res_para_Big_agg_sd_long, digits = c(0,0,0,0,3,3))
print(latex_table, include.rownames = FALSE)

# ------ Permutation Test ------
a_res_Big_agg_sd <- apply(a_res[, , , , , 7, ], c(2,5,6), function(x) sd(x, na.rm=TRUE))
dimnames(a_res_Big_agg_sd) <- list(c(50, 100, 200),
                                   c("Fixed: Lagged", "RE: Lagged", "Fixed: Innovation", "RE: Innovation", "Between"),
                                   c("TPR", "FPR"))
a_res_Big_agg_sd_long <- melt(a_res_Big_agg_sd)
colnames(a_res_Big_agg_sd_long) <- c("#Subjects","ParType", "Outcome", "SD")
n <- 50 # number of repetitions of design
a_res_Big_agg_sd_long$SE <- a_res_Big_agg_sd_long$SD / sqrt(n)
a_res_Big_agg_sd_long$SE <- round(a_res_Big_agg_sd_long$SE, 3)
a_res_Big_agg_sd_long$SD <- round(a_res_Big_agg_sd_long$SD, 3)

# Print the LaTeX table
latex_table <- xtable(a_res_Big_agg_sd_long, digits = c(0,0,0,0,3,3))
print(latex_table, include.rownames = FALSE)

# Across everything for TPR/FPR
round(apply(a_res, 7, function(x) sd(x, na.rm=TRUE)) / sqrt(50), 3) # n=50
round(apply(a_res, 7, function(x) sd(x, na.rm=TRUE)) / sqrt(100), 3) # n=100


# --------------------------------------------------------
# ------- Count Resamples --------------------------------
# --------------------------------------------------------

# VAR lag matrix
resamp1_agg_mean <- apply(a_resamples1, 4, mean)
resamp1_agg_max <- apply(a_resamples1, 4, max)
resamp1_long_mean <- melt(resamp1_agg_mean)
resamp1_long_max <- melt(resamp1_agg_max)
resamp1_long_cmb <- cbind(resamp1_long_mean, resamp1_long_max)

# Residual cov matrix
resamp2_agg_mean <- apply(a_resamples2, 4, mean)
resamp2_agg_max <- apply(a_resamples2, 4, max)
resamp2_long_mean <- melt(resamp2_agg_mean)
resamp2_long_max <- melt(resamp2_agg_max)
resamp2_long_cmb <- cbind(resamp2_long_mean, resamp2_long_max)

# Combine
m1 <- cbind(rep("VAR", 3), v_Delta, resamp1_long_cmb)
colnames(m1) <- c("Matrix", "Delta", "Mean(Resamp)", "(Max Resamp)")
m2 <- cbind(rep("Innovation", 3), v_Delta, resamp2_long_cmb)
colnames(m2) <- colnames(m1)
resamp_all <- rbind(m1, m2)

xtable_output <- xtable(resamp_all, include.rownames = FALSE, digits=c(0,0, 2, 0, 0))
# Print the xtable without row names
print(xtable_output, include.rownames = FALSE)










