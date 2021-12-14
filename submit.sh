#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1

##SBATCH --exclude=markstein.epp.warwick.ac.uk 


# module purge

# module load GCC/8.3.0  OpenMPI/3.1.4 SciPy-bundle/2019.10-Python-3.7.4 

# module load matplotlib/3.1.1-Python-3.7.4

# module load ASE/3.19.1-Python-3.7.4



srun python hans -i 2e4 -w 10 -l 10 -b 1 -c 32 -t 3600 -p 8 >> hsa_ns_output_new_parallel.txt
 

