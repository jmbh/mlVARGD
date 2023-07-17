# jonashaslbeck@protonmail.com; July 17th, 2023

# --------------------------------------------------------
# ------- What is happening here? ------------------------
# --------------------------------------------------------

# Script to reproduce the tutorial shown in the paper

# --------------------------------------------------------
# ------- Load Packages ----------------------------------
# --------------------------------------------------------

# library(devtools)
# install_github("jmbh/mnet")
library(mnet)

library(readr)
library(plyr)

library(mlVAR)

library(qgraph)

# --------------------------------------------------------
# ------- Code -------------------------------------------
# --------------------------------------------------------

l_outlist <- list()

head(dataKoval13)

# --- Run Test ---
output <- mlVAR_GC(data = dataKoval13,
                   vars = colnames(dataKoval13)[6:12],
                   idvar = "PID",
                   dayvar = "UNIT",
                   beepvar = "OCCASION",
                   groups = "group",
                   nP = 1000,
                   contemporaneous = "orthogonal",
                   temporal = "orthogonal",
                   nCores = 32,
                   pbar = TRUE)

# Save

l_outlist$GC <- output

# --- Run Individual Models ---

out_g1 <- mlVAR(data = dataKoval13[dataKoval13$group==1,],
                vars = colnames(dataKoval13)[6:12],
                idvar = "PID",
                dayvar = "UNIT",
                beepvar = "OCCASION",
                contemporaneous = "orthogonal",
                temporal = "orthogonal")

l_outlist$g1 <- out_g1


out_g2 <- mlVAR(data = dataKoval13[dataKoval13$group==2,],
                vars = colnames(dataKoval13)[6:12],
                idvar = "PID",
                dayvar = "UNIT",
                beepvar = "OCCASION",
                contemporaneous = "orthogonal",
                temporal = "orthogonal")

l_outlist$g2 <- out_g2

saveRDS(l_outlist, file="mnet_output_Koval2013.RDS")
l_outlist <- readRDS("Tutorial/mnet_output_Koval2013.RDS")

# --- Check Results ---

output <- l_outlist$GC
out_g1 <- l_outlist$g1
out_g2 <- l_outlist$g2

# output$Runtime_min
# output$Pval$Lagged_fixed
#
# mean(output$Pval$Lagged_fixed<=0.05)
# sum(output$Pval$Lagged_fixed<=0.05)
#
# output$EmpDiffs

# For Tutorial
l_nets <- list()
l_nets$phi_group0 <- out_g1$results$Beta$mean[, , 1]
l_nets$phi_group1 <- out_g2$results$Beta$mean[, , 1]
l_nets$phi_diff <- output$EmpDiffs$Lagged_fixed[, , 1]
l_nets$phi_diffs_sig <- l_nets$phi_diff
l_nets$phi_diffs_sig[output$Pval$Lagged_fixed>0.05] <- 0

output$Pval$Lagged_fixed

titles <- c("CESD Low (<= 16)", "CESD High (> 16)",
            "Differences: Group Low - Group High",
            "Significant Differences")

# --- Make Figure ---

pdf("Figures/Fig_Tutorial_Kovel2013.pdf", width = 7, height = 7)
par(mfrow=c(2,2))
for(i in 1:4) qgraph(t(l_nets[[i]]), # needs transpose, because qgraph plots X[2,1] as 2->1
                     layout = "circle",
                     edge.labels = TRUE,
                     title = titles[i],
                     theme = "colorblind",
                     maximum=0.2,
                     mar=rep(5,4),
                     vsize=13, esize=15, asize=9, edge.label.cex=1.5)

dev.off()
