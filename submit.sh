#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1

##SBATCH --exclude=markstein.epp.warwick.ac.uk 


# module purge

# module load GCC/8.3.0  OpenMPI/3.1.4 SciPy-bundle/2019.10-Python-3.7.4 

# module load matplotlib/3.1.1-Python-3.7.4

# module load ASE/3.19.1-Python-3.7.4



srun python hans -i 2e4 -R -f NS_32_1mer.100.100.1/ -t 28800 -p 1 >> hsa_ns_output.txt
 

