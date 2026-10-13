#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --job-name=composite
#SBATCH --partition=hm
#SBATCH --exclude=rmhm068
#SBATCH --error=composite.%J.err
#SBATCH --output=composite.%J.out
#SBATCH --mem=800

module load openmpi/4.1.4
module load gcc/11.2.0
export OMP_NUM_THREADS=1

cd $SLURM_SUBMIT_DIR

time mpirun -np $SLURM_NTASKS /home/IITB/multiscale-mechanics/ashwani12/software/builds/manybody-LAMMPS/lammps-29Aug2024/src/lmp_mpi -in run_deform_final.in > out.dat
