# Reproducibility Archive

Reproducibility Archive for Paper "Testing for Group Differences in Multi-level Vector Autoregressive Models"; Preprint XXXX


### Main Simulation Study

The files to reproduce the main simulation study reported in the paper can be found in the folder `Simulation_Main`. The files are:

- `aux_Sim.R` contains a script to generate true models and data that is being used in the main simulation script
- `SimMain_Script.R` contains the main simulation script which runs 1 repetition of the design
- `submit_jobs.sh` and `submit_all.sh` are bash scripts that were used to run `SimMain_Script.R` 50 times parallelized over 32 cores on the Snellius supercomputer cluster of the University of Amsterdam; each call of `SimMain_Script.R` produces an output file which can be found in `output` in the `Simulation_Main` folder
- `Evaluation.R` takes the 50 output files from the simulation and creates the three results figures in the main text

The running time per node parallelized over 32 cores was roughly 9 hours.

NOTE: Github currently throws an error when uploading the quite large simulation output files, so they are not actually available yet. I'm working on sorting out this issue.

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









