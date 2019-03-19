#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.

#SBATCH --job-name=snakemake
#SBATCH --output=/home/bjf79/slurm_stdout_logs/$(date +%Y%m%d_%H-%M)_%j.txt
#SBATCH --time=4-00:00:00
#SBATCH --partition=gilad
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

echo "hello world3"

#activate environment
source activate my_ChimpEQTL_env

snakemake --configfile config_test.yaml dedup/test-lane1.sorted.bam dedup/test2-lane1.sorted.bam  --jobs 500 -p --ri --use-conda --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"
