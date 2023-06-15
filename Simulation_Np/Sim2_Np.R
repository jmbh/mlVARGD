# jonashaslbeck@protonmail.com; June 15th, 2023

# --------------------------------------------------------------
# ---------- Get Iteration Number ------------------------------
# --------------------------------------------------------------

# !/usr/bin/env Rscript
iter <- commandArgs(trailingOnly=TRUE)
print(iter)
iter <- as.numeric(iter)

# ----------------------------------------------------------------------
# ----- Load Packages --------------------------------------------------
# ----------------------------------------------------------------------

# Panic Model Sim
library(mlVAR)
library(mixnet)

# Parallelization
library(foreach)
library(parallel)
library(doParallel)


# ----------------------------------------------------------------------
# ----- Load Data --------------------------------------------------
# ----------------------------------------------------------------------

l_data <- readRDS("mlVARGD_sim_p6_Final.RDS")

data_cmb <- rbind(l_data[[1]], l_data[[2]])
data_cmb <- cbind(data_cmb, c(rep(1, nrow(l_data[[1]])),
                              rep(2, nrow(l_data[[2]]))))
colnames(data_cmb)[8] <- "group"

head(data_cmb)

# ----------------------------------------------------------------------
# ----- Iterate --------------------------------------------------------
# ----------------------------------------------------------------------

set.seed(iter)

# ----- Iterating -----
timer_total <- proc.time()[3]

timer <- proc.time()[3]
out_mlVAR_GC <- mlVAR_GC(data = data_cmb,
                         vars = paste0("V", 1:6),
                         idvar = "id",
                         groupo = "group",
                         nP = 10000,
                         verbose = FALSE,
                         saveModels = FALSE,
                         contemporaneous = "orthogonal",
                         temporal = "orthogonal",
                         nCores = 32)
final_time <- proc.time()[3] - timer


# print total time of nodes
print(paste0("Full Timing Iteration ", iter, ":"))
proc.time()[3] - timer_total


# ----------------------------------------------------------------------
# ----- Export ---------------------------------------------------------
# ----------------------------------------------------------------------

# Save
saveRDS(out_mlVAR_GC, file = paste0("mlVARGD_Sim1_N10k_Iter", iter,".RDS"))





