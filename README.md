# hans


[![DOI](https://sandbox.zenodo.org/badge/413444999.svg)](https://sandbox.zenodo.org/doi/10.5072/zenodo.276)


Implementation of a nested sampling algorithm taking advantage of a slightly modified version of David Quigley's hs_alkane package.

One can run a nested sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. Inputs for the program can be fed to it via `stdin`. An example of calling the program may look like:

`python mpihans.py < input.txt > log.out `

An example of what an input file may look like is written below:

<pre>
move_ratio = 3,48,32,0,3,3
nbeads = 6
nchains = 16
nwalkers = 20
iterations = 800000
walklength = 40
time = 190000
analyse = 1
</pre>
Inputs can be commented out for convenience if not wanted, but note that if loading a restart, the quantities in the restart file will overwrite those supplied in the input file.

Running the program will generate a folder which contains: a `restart.hdf5` file, allowing simulations to be restarted from from previous runs, an `volumes.txt` file which is compatible for use with `ns_analyse` from the `pymatnest` package, and a `traj.extxyz` file, allowing one to view the trajectory with the largest volume after every 1000 iterations. 

As the system is athermal, the output isn't truly an energy, but a volume, and therefore it may be more accurate when performing analysis to examine the packing fraction of the chains. For this purpose, the intersecting spheres notebook has been written, which demonstrates how to obtain the volume for a chain, as well as an equation that one can use to calculate the volume for a chain of length N. 

NesSa is a module which also grants the ability to construct slightly different calculations using python, as well as provide some other useful functionalities which a user may require.

## List of input arguments

`analyse` int. Produce a compressibility vs pressure plot of the system once the simulation is finished. Should be 0 or 1.

`bondangle` float. The angle formed by three consecutive spheres within a chain.

`bondlength`  float. The distance between bonds within a chain.

`directory` string. The folder to create if a new run is being started, or the folder to search inside for the restart file if a run is being continued.

`initial_config` string. File to import for starting configurations. This configuration will be cloned and sent to all walkers, then undergoing a brief Monte Carlo walk before the run starts in order to randomise them. Useful if starting from particular structures such as ringed alkanes.

`min_aspect_ratio` float. Smallest allowed distance between parallel faces for cell normalised to unit volume. A higher value restricts the system to more cube-like cell shapes. Should be between 0 and 1.

`move_ratio` 6 floats separated by commas. Ratio of moves to use when performing Monte Carlo walks. Values correspond with "volume moves", "translational moves", "rotational moves", "dihedral moves", "shear moves", "stretch moves".

`nbeads` int. The number of beads per chain.

`nchains` int. The number of the chains in each simulation cell

`restart_file` string. The file from which to restart a run from.

`walklength` int. The number of "sweeps" performed per iteration on each cpu, constituting a Monte Carlo walk. A sweep is defined as a number of Monte Carlo moves which should change each degree of freedom within the system once on average.

