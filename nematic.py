import numpy as np
import argparse as ap
import ase.io
import matplotlib.pyplot as plt

def parse_args():

    """Parse command line arguments. Returns an empty Namespace containing the attributes."""

    parser = ap.ArgumentParser()
    parser.add_argument("atoms_file",metavar = "f",type=str, help =
                        "Which config file to calculate nematic order param for.\n")
    parser.add_argument("-b","--beads",type=int, help =
                        "How many beads per chain in the configurations.\n")
    parser.add_argument("-c","--chains",type=int, help =
                        "How many chains there are total in the configuration.\n")
    parser.add_argument("-P","--pick_configs",action='store_true', help =
                        "Whether to pick out configs with a nematic order parameter in a certain range")
    parser.add_argument("-l","--lbound",type=float, help =
                        "Lower bound for nematic order parameter for configs to be picked.\n")
    parser.add_argument("-u","--ubound",type=float, help =
                        "Upper bound for nematic order parameter for configs to be picked.\n")
                        
    return parser.parse_args()


def nematic_order_param(configs,Nbeads,Nchains):
    """
    
    Plots the nematic order parameter over time for an array of Atoms Objects or a trajectory.
    
    configs: Array of Atoms Objects/Trajectory for the atoms objects which needs to be calculated
    Nbeads: Number of beads in the chain which is being calculated
    Nchains: The number of chains present in the system which is being calculated
    
    """

    order_param_array=[]
    Nconfigs = len(configs) 


    for img in configs:
        
        Q_sum = np.zeros((3,3)) # Initialise Q_sum
        
        for i in range(Nchains):
            
            # Calculate the vector from the beginning chain to the end chain
            # -- Note that I'm letting ASE just do all the work for me here!

            a0 = i*Nbeads      # index of zeroth bead on chain i
            a1 = i*Nbeads + Nbeads-1  # index of last bead on chain i 

            # Find end-to-end vector
            chainvector = img.get_distance(a0, a1, mic=True, vector=True)
            
            # Alternatively use molecular direction from inertia tensor
            #chainvector = molecular_axis(img, Nbeads, i)
            
            # Normalise chainvector
            chainvector = chainvector/np.linalg.norm(chainvector)

            # Accumulate Q_sum
            Q_sum = Q_sum + 1.5 * np.outer(chainvector, chainvector) - 0.5 * np.identity(3)
            
           
        # Find the average Q
        Qav = Q_sum/Nchains
        
        # Diagonalise this
        eigenvals, eigenvects = np.linalg.eig(Qav)
        
        # The maximum positive eigenvalue is the P2 order parameter
        P2 = np.amax(eigenvals)
        
        order_param_array.append(P2)
        
        # If we wanted the director then this is the corresponding eigenvector
        #k = np.where(eigenvals == P2)
        #director = eigenvects[k]
    return np.array(order_param_array)

def plot_fig(array):
    x_vals = range(len(array))
    plt.figure(figsize = (8,4),dpi = 100)
    plt.plot(x_vals,array, lw=0.4)
    plt.scatter(x_vals,array, marker = 'x')
    plt.xlabel("Config no.")
    plt.ylabel("Nematic Order Parameter")
    plt.tight_layout()
    plt.savefig("nematic.png")
    return

def config_picker(atoms,nem_array,lower_bound,upper_bound):
    
    picked_indices=np.where((nem_array>lower_bound)&(nem_array<upper_bound))
    picked_atoms = [atoms[i] for i in picked_indices[0]]
    ase.io.write("picked_nematic.extxyz", picked_atoms)
    return


def nematic_analysis():
    args = parse_args()
    configs = ase.io.read(args.atoms_file,index=":")
    nematic_array = nematic_order_param(configs,args.beads,args.chains)
    if args.pick_configs and (args.lbound is None or args.ubound is None):
        args.error("--pick_configs requires --lbound and --ubound.")
    plot_fig(nematic_array)
    for i,nem in enumerate(nematic_array):
        print(i,nem)
    if args.pick_configs:
        config_picker(configs,nematic_array,args.lbound,args.ubound)
    return

nematic_analysis()