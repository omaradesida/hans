#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2012
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=hetmathsys
##SBATCH --signal=B:TERM@30
 


# module purge

# module load GCC/8.3.0  OpenMPI/3.1.4 SciPy-bundle/2019.10-Python-3.7.4 

# module load matplotlib/3.1.1-Python-3.7.4

# module load ASE/3.19.1-Python-3.7.4



srun -n 4 python mpihans.py<input.txt #> MPI_hans.out
# srun --signal=B:USR1@30 python signal_test.py > signal_test.out







