#!/bin/bash
#SBATCH --job-name=ssages_test
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --time=00:15:00
#SBATCH --partition=boltzmann
#SBATCH --qos=debug
#SBATCH --ntasks=12
#SBATCH --share
#SBATCH --mail-user=bezik@uchicago.edu
#SBATCH --mail-type=ALL

mpiexec -np 12 ./ssages FTS.json

