#!/bin/bash
#SBATCH --job-name ob_prior
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7000m
#SBATCH --time=00:15:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard
#SBATCH --mail-type=NONE
#SBATCH --export=ALL
#SBATCH --output=%x-%j.log
#SBATCH --array=1-1000

cd $SLURM_SUBMIT_DIR
Rscript --verbose poisson_sample_size_log_theta.R $SLURM_ARRAY_TASK_ID 0.65
