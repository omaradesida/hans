import numpy as np
import argparse as ap
import ase.io
import matplotlib.pyplot as plt

def parse_args():

    """Parse command line arguments. Returns an empty Namespace containing the attributes."""

    parser = ap.ArgumentParser()
    parser.add_argument("atoms_file",metavar = "f",type=str, help =
                        "Which config file to calculate order param for.\n")
    parser.add_argument("-m","--mode",type=str, help =
                        "Which order param to calculate:\n \
                        'nematic'('n') or 'translational'('t').\n", default = "nematic")
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


def trans_order_param(configs,Nbeads,Nchains):
    """
    configs: array of configs to calculate order parameter for. Can be generated using ase.io.read
    Nbeads : Number of beads per chain in the configuration
    Nchains: Number of chains in the configuration"""
    order_param_array = []
    
    for config in configs[:]:
        rho_sum = []     

        recip_cell = 2*np.pi*config.cell.reciprocal()
        #com = centre_of_masses(config,Nbeads,Nchains, fractional = True)
        r = centre_of_masses(config,Nbeads,Nchains, fractional = True)

        l = (config.get_volume()**(1/3))
        k = recip_cell
        # k = recip_cell/l
        for i in range(Nchains*Nbeads):
            rho = np.exp(1j*2*np.pi*np.sum(r[i]))
            
            rho_sum.append(rho)
        rho_sum = np.array(rho_sum)
        rho_sum = np.sum(rho_sum,axis=0)/(Nchains)
        param = np.abs(rho_sum)
        order_param_array.append(param)
    return np.array(order_param_array)

def centre_of_masses(config,Nbeads,Nchains, fractional = True):
    """
    Calculate the centre of mass for a chain in the configuration.
    """

        
    com_array = []
    if fractional:
        for i in range(Nchains):
            a0 = i*Nbeads      # index of zeroth bead on chain i
            a1 = i*Nbeads + Nbeads-1  # index of last bead on chain i
            if Nbeads==1:
                com_array.append(config.get_scaled_positions()[a0])
            else:
                com_array.append(config.get_scaled_positions()[a0:a1].mean(0))
    else:
        for i in range(Nchains):
            a0 = i*Nbeads      # index of zeroth bead on chain i
            a1 = i*Nbeads + Nbeads-1  # index of last bead on chain
            if Nbeads==1:
                com_array.append(config.get_positions()[a0])
            else:
                com_array.append(config.get_positions()[a0:a1].mean(0))

    return com_array

def nematic_order_param(configs,Nbeads,Nchains):
    """
    
    Calculate the nematic order parameter for each config for an array of Atoms Objects or a trajectory.
    
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
    try:
        plt.scatter(x_vals,array, marker = 'x')
    except:
        plt.scatter(x_vals,array[:,0], marker = 'x')
        plt.scatter(x_vals,array[:,1], marker = 'x')
        plt.scatter(x_vals,array[:,2], marker = 'x')


    plt.xlabel("Config no.")
    plt.ylabel("Order Parameter")
    plt.tight_layout()
    plt.savefig("order_param.png")
    return

def config_picker(atoms,nem_array,lower_bound,upper_bound):
    
    picked_indices=np.where((nem_array>lower_bound)&(nem_array<upper_bound))
    picked_atoms = [atoms[i] for i in picked_indices[0]]
    ase.io.write("picked_order_param.extxyz", picked_atoms)
    return


def order_param_analysis():
    args = parse_args()
    configs = ase.io.read(args.atoms_file,index=":")
    if args.mode in ["nematic","n"]:
        order_param = nematic_order_param
    elif args.mode in ["translational","t"]:
        order_param = trans_order_param

    order_param_array = order_param(configs,args.beads,args.chains)
    if args.pick_configs and (args.lbound is None or args.ubound is None):
        args.error("--pick_configs requires --lbound and --ubound.")
    for i,nem in enumerate(order_param_array):
        print(i,nem)
    plot_fig(order_param_array)
    if args.pick_configs:
        config_picker(configs,order_param_array,args.lbound,args.ubound)
    return

if __name__ == "__main__":
    order_param_analysis()