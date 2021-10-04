# hsa_NS
Implementation of a nested sampling algorithm taking advantage of David Quigley's hs_alkane package.

One can run a nested_sampling calculation of chains using hard sphere potentials, using the `hs_alkane` package to perform Monte Carlo (MC) moves as well as
store the particle positions. One can run a nested sampling calculation using the following command:

`python hs_alkane_NS.py N_walkers N_chains N_beads N_iter steps_per_walk`

using appropriate values for the arguments.

`N_walkers`: How many walkers should be used for each run

`N_chains`: How many chains should be used per simulation box.

`N_beads`: How many beads there are in each chain. Note that a chain of length one will produce a box populated with spheres.

`N_iter`: How many iterations of the nested sampling algorithm to be ran

`steps_per_walk`: How many steps are performed when a random walker is cloned. Steps use MC moves to generate new configurations.
