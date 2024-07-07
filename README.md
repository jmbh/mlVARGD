# Reproducibility Archive

Reproducibility Archive for Paper "Testing for Group Differences in Multi-level Vector Autoregressive Models"; Preprint https://psyarxiv.com/dhp8s/


### Main Simulation Study

The files to reproduce the main simulation study reported in the paper can be found in the folder `Simulation_Main`. The files are:

- `aux_Sim.R` contains a script to generate true models and data that is being used in the main simulation script
- `SimMain_Script.R` contains the main simulation script which runs 1 repetition of the design
- `submit_jobs.sh` and `submit_all.sh` are bash scripts that were used to run `SimMain_Script.R` 50 times parallelized over 32 cores on the Snellius supercomputer cluster of the University of Amsterdam; each call of `SimMain_Script.R` produces an output file which can be found in `output` in the `Simulation_Main` folder
- `Evaluation.R` takes the 50 output files from the simulation and creates the three results figures in the main text

The running time per node parallelized over 32 cores was roughly 9 hours.

NOTE: The simulation results are 2.28GB, which is too large to share on Github. However, I am happy to share the files upon request.

### Simulation Study on Np in Appendix

Appendix A reports a small simulation study to illustrate the number of permutations needed to approximate different p-values. The files to reproduce this simulation are in the folder ``Simulation_Np``.

- `Sim2_Np.R` contains the script to run the permutation test 10 times with 10,000 repetitions on an example data set with six variables
- `mlVARGD_sim_p6_Final.RDS` is the synthetic example dataset used
- `submit_jobs.sh` and `submit_all.sh` are bash scripts used as in the main simulation
- `Evaluation.R` takes the 10 output file and creates the figure in the appendix


### Tutorial

The folder `Tutorial` contains files to reproduce the analysis and results in the tutorial section.

- `Tutorial_Koval13.R` contains the script for the tutorial, parts of which is also shown in the paper. It also contains the code to produce the figure shown in the tutorial section.
- `mnet_output_Koval2013.RDS` is the output file from the permutation test run in the tutorial section. We provide it separately here, because the script runs relatively long



### Session Info

sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.6 (Ootpa)

Matrix products: default
BLAS/LAPACK: /gpfs/admin/_hpc/sw/arch/AMD-ZEN2/Centos8/EB_production/2021/software/FlexiBLAS/3.0.4-GCC-10.3.0/lib64/libflexiblas.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] qgraph_1.9.8       plyr_1.8.9         readr_1.4.0        doParallel_1.0.17 
 [5] iterators_1.0.14   foreach_1.5.2      mlVAR_0.5.3        corpcor_1.6.10    
 [9] mnet_0.1.3         devtools_2.4.1     usethis_2.0.1      xtable_1.8-4      
[13] reshape2_1.4.4     RColorBrewer_1.1-3 scales_1.3.0      

loaded via a namespace (and not attached):
 [1] minqa_1.2.7             colorspace_2.1-0        ellipsis_0.3.2         
 [4] htmlTable_2.4.2         base64enc_0.1-3         fs_1.6.4               
 [7] rstudioapi_0.16.0       lavaan_0.6-18           remotes_2.4.0          
[10] fansi_1.0.6             mvtnorm_1.2-5           codetools_0.2-18       
[13] splines_4.1.0           mnormt_2.1.1            cachem_1.1.0           
[16] knitr_1.47              glasso_1.11             texreg_1.39.3          
[19] pkgload_1.3.4           Formula_1.2-5           nloptr_2.1.0           
[22] cluster_2.1.2           png_0.1-8               compiler_4.1.0         
[25] httr_1.4.7              backports_1.5.0         Matrix_1.3-4           
[28] fastmap_1.2.0           cli_3.6.3               htmltools_0.5.8.1      
[31] tools_4.1.0             igraph_2.0.3            coda_0.19-4.1          
[34] gtable_0.3.5            glue_1.7.0              clusterGeneration_1.3.8
[37] dplyr_1.1.4             Rcpp_1.0.12             vctrs_0.6.5            
[40] nlme_3.1-152            psych_2.4.3             xfun_0.45              
[43] stringr_1.5.1           proto_1.0.0             lme4_1.1-35.2          
[46] lifecycle_1.0.4         gtools_3.9.5            MASS_7.3-54            
[49] graphicalVAR_0.3.4      hms_1.1.3               memoise_2.0.1          
[52] pbapply_1.7-2           gridExtra_2.3           ggplot2_3.5.1          
[55] pander_0.6.5            rpart_4.1-15            stringi_1.8.4          
[58] fastDummies_1.7.3       checkmate_2.3.1         boot_1.3-28            
[61] pkgbuild_1.4.4          shape_1.4.6.1           rlang_1.1.4            
[64] pkgconfig_2.0.3         arm_1.14-4              evaluate_0.24.0        
[67] lattice_0.20-44         purrr_0.3.4             htmlwidgets_1.6.4      
[70] tidyselect_1.2.1        magrittr_2.0.3          R6_2.5.1               
[73] generics_0.1.3          Hmisc_5.1-3             gsubfn_0.7             
[76] pillar_1.9.0            foreign_0.8-81          withr_3.0.0            
[79] survival_3.2-11         abind_1.4-5             nnet_7.3-16            
[82] tibble_3.2.1            fdrtool_1.2.17          utf8_1.2.4             
[85] MplusAutomation_1.1.1   rmarkdown_2.27          jpeg_0.1-10            
[88] grid_4.1.0              data.table_1.15.4       pbivnorm_0.6.0         
[91] digest_0.6.35           stats4_4.1.0            munsell_0.5.1          
[94] glmnet_4.1-8            sessioninfo_1.1.1       quadprog_1.5-8    








