# hsa_NS
Implementation of a nested sampling algorithm taking advantage of a slightly modified version of David Quigley's hs_alkane package.

One can run a nested sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. One can run a nested sampling calculation using the following command:

`python hs_alkane_NS.py [options]`

using appropriate values for the arguments.
<pre>  -h, --help            show this help message and exit
  -i --iterations,      ITERATIONS
                        Number of iterations to run
  -w --nwalkers,        NWALKERS
                        Number of walkers to use
  -c --nchains,         NCHAINS
                        Number of chains in the simulation box
  -b --nbeads,          NBEADS
                        Number of beads in each chain
  -l --walklength,      WALKLENGTH
                        Number of sweep for each random walk step
  -R, --restart,        Whether or not to restart from a previous attempt
  
  -f --restart_folder,  RESTART_FOLDER
                        Configurations file to restart from. Energies and Trajectory file should have the same file naming convention. 
                        i.e if restart file is &apos;foo.restart&apos; energies and trajectory file should be foo.energies and
                        foo.extxyz
</pre>

Simulations will generate an energies file which is compatible for use with `ns_analyse` from the `pymatnest` package. Note that as the system is athermal, the output isn't truly an energy and therefore it may be more accurate when performing analysis to examine the packing fraction of the chains. For this purpose, the intersecting spheres notebook has been written, which demonstrates how to obtain the volume for a chain, as well as an equation that one can use to calculate the volume for a chain of length $N$.
