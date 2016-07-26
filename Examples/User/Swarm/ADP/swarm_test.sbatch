#!/bin/bash
#SBATCH --job-name=ADP_test
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=boltzmann
#SBATCH --ntasks=22
#SBATCH --share
#SBATCH --mail-user=bezik@uchicago.edu
#SBATCH --mail-type=ALL

mpiexec -np 22 ./ssages Swarm.json

