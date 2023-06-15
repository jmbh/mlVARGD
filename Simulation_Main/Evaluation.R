# jonashaslbeck@protonmail.com; June 15th, 2023

# --------------------------------------------------------
# ------- What is happening here? ------------------------
# --------------------------------------------------------

# Analyze output of Simulation Study 1 on Permutations


# --------------------------------------------------------
# ------- Load Packages ----------------------------------
# --------------------------------------------------------

library(scales)
library(RColorBrewer)

# --------------------------------------------------------
# ------- Some Aux for Evaluation ------------------------
# --------------------------------------------------------

plotLabel <- function(x, srt=0, col="black",
                      xpos=.6, ypos=.6, cex=1.4) {
  par(mar=rep(0, 4))
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  text(xpos, ypos, x, srt=srt, cex=cex, col=col)
}

# --------------------------------------------------------
# ------- Loading Sim Results ----------------------------
# --------------------------------------------------------

results_folder <- "output/"

l_all <- list()
v_files <- list.files(paste0("Simulation_Main/", results_folder))
n_files <- length(v_files)
l_files <- list()
for(i in 1:n_files) l_files[[i]] <- readRDS(paste0("Simulation_Main/", results_folder, v_files[i]))


date_of_analysis <- "June15th_23"

# --------------------------------------------------------
# ------- Some Sanity Checks -----------------------------
# --------------------------------------------------------

# --------------------------------------------------------
# ------- Preprocess Results: Permutation Test -----------
# --------------------------------------------------------

nIter <- n_files
n_delta <- 3
n_Nsubj <- 3
n_prop <- 2

n_pvals <- 100
v_pvals <- seq(0, 1, length=n_pvals)

## Storage
# iterations, n_subj, n_prop, n_delta, 5 vartypes, 100 x pval thresh, 2 outcomes (TPR, FPR)
a_res <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta, 5, n_pvals, 2))

for(i in 1:nIter) {
  j <- 1   # counter for design matrix
  for(i_sub in 1:3) { # number of subjects
    for(i_pr in 1:2) { # proportion
      for(i_d in 1:3) { # delta

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

        j <- j+1 # update counter; should be 18+1 in the end
      }
    }
  }
}

# Aggregate across repetitions of simulation study
a_res_agg <- apply(a_res, 2:7, mean)

dim(a_res)


# --------------------------------------------------------
# ------- Main Results Figure: ROC Curves Permutation ----
# --------------------------------------------------------

library(RColorBrewer)
cols <- brewer.pal(3, "Set1")


pdf(paste0("Figures/Fig_Sim1_",date_of_analysis,".pdf"), width=7, height=9)

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
plotLabel(expression(paste(Delta, " = 0.05")))
plotLabel(expression(paste(Delta, " = 0.10")))
plotLabel(expression(paste(Delta, " = 0.25")))
# Left Col
plotLabel("Lagged Fixed", srt=90)
plotLabel("Lagged Var", srt=90)
plotLabel("Resid Fixed", srt=90)
plotLabel("Resid Var", srt=90)
plotLabel("Between", srt=90)

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
        for(alpha in 1:3) points(a_res_agg[i_sub, i_pr, i_d, type, c(2,6,11)[alpha], 2],
                                 a_res_agg[i_sub, i_pr, i_d, type, c(2,6,11)[alpha], 1],
                                 col = cols[i_sub],
                                 pch=alpha)
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


# --------------------------------------------------------
# ------- Preprocess Results: Parametric Test ------------
# --------------------------------------------------------

nIter <- n_files
n_delta <- 3
n_Nsubj <- 3
n_prop <- 2

n_pvals <- 100
v_pvals <- seq(0, 1, length=n_pvals)

## Storage
# iterations, n_subj, n_prop, n_delta, 3 para types, 100 x pval thresh, 2 outcomes (TPR, FPR)
a_res_para <- array(NA, dim=c(nIter, n_Nsubj, n_prop, n_delta, 3, n_pvals, 2))

for(i in 1:nIter) {
  j <- 1   # counter for design matrix
  for(i_sub in 1:3) {
    for(i_pr in 1:2) {
      for(i_d in 1:3) {

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

        j <- j+1 # update counter; should be 18+1 in the end
      }
    }
  }
}

# Aggregate across repetitions of simulation study
a_res_para_agg <- apply(a_res_para, 2:7, mean)

dim(a_res_para)


# --------------------------------------------------------
# ------- Main Results Figure: ROC Curves Parametric -----
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

v_pvals <- seq(0, 1, length=n_pvals)


pdf(paste0("Figures/Fig_Sim1_Parametric_",date_of_analysis,".pdf"), width=7, height=9*(3/5))

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
plotLabel(expression(paste(Delta, " = 0.05")))
plotLabel(expression(paste(Delta, " = 0.10")))
plotLabel(expression(paste(Delta, " = 0.25")))
# Left Col
plotLabel("Lagged Fixed", srt=90)
plotLabel("Resid Fixed", srt=90)
plotLabel("Between", srt=90)

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
        for(alpha in 1:3) points(a_res_para_agg[i_sub, i_pr, i_d, type, c(2,6,11)[alpha], 2],
                                 a_res_para_agg[i_sub, i_pr, i_d, type, c(2,6,11)[alpha], 1],
                                 col = cols[i_sub],
                                 pch=alpha)

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


# --------------------------------------------------------
# ------- Figure: Difference Between Para vs. Non-Para ---
# --------------------------------------------------------

cols <- brewer.pal(3, "Set1")

v_pvals <- seq(0, 1, length=n_pvals)

# Get subset of permutation results
a_res_nonpara_subset <- a_res_agg[, , , c(1,3,5), , ]

pdf(paste0("Figures/Fig_Sim1_Para_vs_NonPara_",date_of_analysis,".pdf"), width=7, height=9*(3/5))

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
plotLabel(expression(paste(Delta, " = 0.05")))
plotLabel(expression(paste(Delta, " = 0.10")))
plotLabel(expression(paste(Delta, " = 0.25")))
# Left Col
plotLabel("Lagged Fixed", srt=90)
plotLabel("Resid Fixed", srt=90)
plotLabel("Between", srt=90)

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

      if(i_pr==1) a_res_disp <- a_res_nonpara_subset # results permutation test
      if(i_pr==2) a_res_disp <- a_res_para_agg # results parametric test

      # Average across balanced & unbalanced
      a_res_disp_aob <- apply(a_res_disp, c(1,3:6), mean)

      # balance <- 2

      for(i_sub in 1:3) {

        lines(a_res_disp_aob[i_sub, i_d, type, , 2], a_res_disp_aob[i_sub, i_d, type, , 1],
              col=scales::alpha(cols[i_sub], alpha=0.5), lty=c(1,4)[i_pr], lwd=2)


        # Add alphas for 0.01, 0.05, 0.10
        for(alpha in 1:3) points(a_res_disp_aob[i_sub, i_d, type, c(2,6,11)[alpha], 2],
                                 a_res_disp_aob[i_sub, i_d, type, c(2,6,11)[alpha], 1],
                                 col = cols[i_sub],
                                 pch=alpha)
      }
    }

    # Legends
    if(type==1 & i_d==1) {
      cex_size <- 0.9
      legend(.53, .65, legend=paste0("Nsubj = ", c(50, 100, 200)), bty="n", text.col=cols, cex=cex_size)
      legend(.4, .3, legend=c("Non-Parametric", "Parametric"), bty="n", lty=c(1,4), cex=0.8)
    }

    if(type==1 & i_d==2) {
      cex_size <- 0.9
      legend(.5, .45, legend=c(expression(paste(alpha==0.01)),
                                expression(paste(alpha==0.05)),
                                expression(paste(alpha==0.10))), bty="n", pch=1:3, cex=cex_size)
    }


  }
}

dev.off()



