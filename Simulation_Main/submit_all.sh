#!/bin/bash

mkdir "$TMPDIR"/mlVARGD_Sim1/

cd "$HOME"/mlVARGD_Sim1

sbatch -a 1-50 submit_jobs.sh



