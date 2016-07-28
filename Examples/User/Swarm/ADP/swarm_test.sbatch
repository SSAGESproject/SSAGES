#!/bin/bash
#SBATCH --job-name=ADP_test
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=boltzmann
#SBATCH --ntasks=12
#SBATCH --share
#SBATCH --mail-user=bezik@uchicago.edu
#SBATCH --mail-type=ALL

mpiexec -np 12 ./ssages Swarm.json

