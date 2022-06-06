# hans
Implementation of a nested sampling algorithm taking advantage of a slightly modified version of David Quigley's hs_alkane package.

One can run a nested sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. Inputs for the program can be fed to it via `stdin`. An example of calling the program may look like

`python mpihans.py < input.txt > log.out `

The program will create an appropriate directory store simulation data, trajectories and restarts.

using appropriate values for the arguments.

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

Running the program will generate a folder which contains: a `restart.hdf5` file, allowing simulations to be restarted from from previous runs, an `volumes.txt` file which is compatible for use with `ns_analyse` from the `pymatnest` package, and a `traj.extxyz` file, allowing one to view the trajectory with the largest volume after every 1000 iterations. 

As the system is athermal, the output isn't truly an energy, but a volume, and therefore it may be more accurate when performing analysis to examine the packing fraction of the chains. For this purpose, the intersecting spheres notebook has been written, which demonstrates how to obtain the volume for a chain, as well as an equation that one can use to calculate the volume for a chain of length N. 

NesSa is a module which also grants the ability to construct slightly different calculations using python, as well as provide some other useful functionalities which a user may require.


