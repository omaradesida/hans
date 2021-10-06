# hsa_NS
Implementation of a nested sampling algorithm taking advantage of David Quigley's hs_alkane package.

One can run a nested_sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. One can run a nested sampling calculation using the following command:

`python hs_alkane_NS.py [options]`

using appropriate values for the arguments.
<pre>  -h, --help            show this help message and exit
  -i ITERATIONS,         --iterations ITERATIONS
                        Number of iterations to run
  -w NWALKERS,          --nwalkers NWALKERS
                        Number of walkers to use
  -c NCHAINS,           --nchains NCHAINS
                        Number of chains in the simulation box
  -b NBEADS,            --nbeads NBEADS
                        Number of beads in each chain
  -l WALKLENGTH,        --walklength WALKLENGTH
                        Number of sweep for each random walk step
  -R, --restart         Whether or not to restart from a previous attempt
  
  -f RESTART_FOLDER,    --restart_folder RESTART_FOLDER
                        Configurations file to restart from. Energies and Trajectory file should have the same file naming convention. 
                        i.e if restart file is &apos;foo.restart&apos; energies and trajectory file should be foo.energies and
                        foo.extxyz
</pre>
