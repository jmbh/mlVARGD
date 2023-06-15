# jonashaslbeck@protonmail.com; June 1st, 2023

# --------------------------------------------------------------
# ---------- Get Iteration Number ------------------------------
# --------------------------------------------------------------

# !/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
print(iter)
iter <- as.numeric(iter)


# ----------------------------------------------------------------------
# ----- Source & Load Packages -----------------------------------------
# ----------------------------------------------------------------------

# library(devtools)
# install_github("jmbh/mnet")

# Source
source("aux_Sim.R")

# Estimation
library(mnet)
library(corpcor)
library(mlVAR)

# Parallelization
library(foreach)
library(parallel)
library(doParallel)

# ----------------------------------------------------------------------
# ----- Simulating one Iteration ---------------------------------------
# ----------------------------------------------------------------------

# -----------------------
# ----- Storage ---------
# -----------------------

l_results_iter <- list("design" = NULL,
                       "datagen" = list(),
                       "outcomes" = list(),
                       "outcomes_para" = list())


# -----------------------
# ----- Define Design ---
# -----------------------

v_tot_sub <- c(50, 100, 200) # identifiability issue solved!
v_balance <- c(0.5, 0.2) # percentage in first group
v_Delta <- c(0.05, 0.15, 0.25)

# Create design grid
design_grid <- expand.grid(v_tot_sub, v_balance, v_Delta)
colnames(design_grid) <- c("Nsubj", "Balance", "Delta")
n_design <- nrow(design_grid)
l_results_iter$design <- design_grid

v_timer_i <- rep(NA, n_design)

# Total Timing
timer_total <- proc.time()[3]

# Loop through Design
for(i in 1:n_design) {

  # -----------------------
  # ----- Generate Data ---
  # -----------------------

  timer_data_i <- proc.time()[3]

  data_gen <- SimData(p = 8, # Fixed, otherwise infeasible
                      Nsubj = design_grid$Nsubj[i],
                      prop = design_grid$Balance[i],
                      Delta = design_grid$Delta[i],
                      Nt = 100,
                      seed = iter, # Keep constant across design, vary across repetitions
                      var_burnin = 50,
                      pbar = FALSE)

  # saveRDS(data_gen, "data_design4.RDS")

  # print total time of nodes
  print(paste0("Timing DataGen Iter: ", i, ":"))
  print(proc.time()[3] - timer_data_i)

  l_results_iter$datagen[[i]] <- data_gen

  # Make grouping variable
  data_cmd <- rbind(data_gen$data[[1]], data_gen$data[[2]])
  groups <- c(rep(1, nrow(data_gen$data[[1]])),
              rep(2, nrow(data_gen$data[[2]])))

  data_cmd <- cbind(data_cmd, groups)

  # ------------------------------------
  # ----- Testing: Permutation ---------
  # ------------------------------------

  timer_est_i <- proc.time()[3]

  out_mlVAR_GC <- mlVAR_GC(data = data_cmd,
                           vars = paste0("V", 1:8), # 8 variables!
                           idvar = "id",
                           groups = "groups",
                           test = "permutation",
                           nP = 1000, # 64; fixed=1000
                           verbose = FALSE,
                           saveModels = FALSE,
                           contemporaneous = "orthogonal",
                           temporal = "orthogonal",
                           nCores = 32,  # 32; runs quicker than 124 for some reason that is unresolved by Surf Support
                           pbar = FALSE)

  l_results_iter$outcomes[[i]] <- out_mlVAR_GC

  # print total time of nodes
  print(paste0("Timing PermTest Iter: ", i, ":"))
  v_timer_i[i] <- proc.time()[3] - timer_est_i
  print(v_timer_i[i])

  # ------------------------------------
  # ----- Testing: Parametric ----------
  # ------------------------------------

  out_mlVAR_GC_para <- mlVAR_GC(data = data_cmd,
                                vars = paste0("V", 1:8), # 8 variables!
                                idvar = "id",
                                groups = "groups",
                                test = "parametric",
                                verbose = FALSE,
                                saveModels = FALSE,
                                contemporaneous = "orthogonal",
                                temporal = "orthogonal",
                                pbar = FALSE)

  l_results_iter$outcomes_para[[i]] <- out_mlVAR_GC_para

  print(i)

} # end loop: through design

l_results_iter$v_timer_i <- v_timer_i

# print Design grid + timing
cbind(design_grid, v_timer_i)


# print total time of nodes
print(paste0("Full Timing Iteration ", iter, ":"))
proc.time()[3] - timer_total

# ----------------------------------------------------------------------
# ----- Export ---------------------------------------------------------
# ----------------------------------------------------------------------

# Save
saveRDS(l_results_iter, file = paste0("mlVARGD_Sim2_Iter", iter,".RDS"))



