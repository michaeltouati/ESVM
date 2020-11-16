#!/bin/bash

#SBATCH --mem=246000		
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16 
#SBATCH -A cin_staff
#SBATCH -p m100_usr_prod
#SBATCH --gres=gpu:4
#SBATCH --time=10:00:00
#SBATCH --job-name="cern_scenarioA_0.18cm"

##SBATCH --mail-type=ALL
##SBATCH --mail-user=n.shukla@cineca.it

#SBATCH --output=job.out.%j
#SBATCH --error=job.err.%j

export OMP_NUM_THREADS=16
export OMP_PLACES=threads

mpirun -n 4 --map-by socket:PE=4 ./esvm input-deck | tee os.stdout

exit










