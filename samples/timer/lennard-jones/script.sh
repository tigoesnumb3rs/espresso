#!/bin/bash
#SBATCH -J lj-script
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tiwa@icp.uni-stuttgart.de

module load fftw/3.3.4/gcc
module load cuda/7.5
module load hdf5/1.8.15/gcc
module load openmpi/1.10.1
module load binutils/2.25
module load gcc/4.9.3
module load boost/1.60.0
module load python/2.7

mpirun /home/tiwa/espresso/build-lj/Espresso lj_liquid.tcl 
