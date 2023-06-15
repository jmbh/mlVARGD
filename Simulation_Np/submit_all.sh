#!/bin/bash

mkdir "$TMPDIR"/mlVARGD_Sim2/

cd "$HOME"/mlVARGD_Sim2

sbatch -a 1-10 submit_jobs.sh



