#!/bin/bash

mkdir "$TMPDIR"/mlVARGD_Sim13/

cd "$HOME"/mlVARGD_Sim13

sbatch -a 1-50 submit_jobs.sh


