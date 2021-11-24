# hsa_NS
Implementation of a nested sampling algorithm taking advantage of a slightly modified version of David Quigley's hs_alkane package.

One can run a nested sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. One can run a nested sampling calculation using the following command:

`python hs_alkane_NS.py [options]`

using appropriate values for the arguments.
<pre>  -h, --help            show this help message and exit
  -i --iterations,      ITERATIONS
                        Number of iterations to run
  -t --time             TIME
                        How much time has been allocated for the program to run
  -w --nwalkers,        NWALKERS
                        Number of walkers to use
  -c --nchains,         NCHAINS
                        Number of chains in the simulation box
  -b --nbeads,          NBEADS
                        Number of beads in each chain
  -l --walklength,      WALKLENGTH
                        Number of sweep for each random walk step
  -p --processes        PROCESSES
                        Number of processes to use in parallel when performing random walks
  -d --bondlength       BONDLENGTH
                        Distance between hard spheres in a chain


  -R, --restart,        Whether or not to restart from a previous attempt
  
  -f --restart_folder,  RESTART_FOLDER
                        Folder where the restart file is located. In the folder, a restart.hdf5 file must be present in order to resume the simulation.
</pre>

Running the program will generate a folder which contains: a `restart.hdf5` file, allowing simulations to be restarted from from previous runs, an `energies.txt` file which is compatible for use with `ns_analyse` from the `pymatnest` package, and a `traj.extxyz` file, allowing one to view the trajectory with the largest volume after every 1000 iterations. Note that the `energies.txt` needs the option `s -1` when using `ns_analyse` as the first line contains metadata about the system, such as the number of beads per chain and the number of walkers. as the system is athermal, the output isn't truly an energy, but a volume, and therefore it may be more accurate when performing analysis to examine the packing fraction of the chains. For this purpose, the intersecting spheres notebook has been written, which demonstrates how to obtain the volume for a chain, as well as an equation that one can use to calculate the volume for a chain of length N. 


