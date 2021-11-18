

from timeit import default_timer as timer


t0 = timer()

from NS_hsa import *
from simulation_parameters import SimulationParameters

#parsing inputs

#parent_dir = os.getcwd()
parent_dir = "." #temporary until hs_alkane can take longer filenames

parameters = SimulationParameters.parse_args()

ns_data = ns_info(parameters)

ns_data.set_directory(parameters.directory)

energies_file = open(ns_data.energies_filename, "a+")

#setting params for MC_run
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']
ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5
move_ratio = np.zeros(6)
move_ratio[ivol] = 1
move_ratio[itrans] = 3.0*ns_data.nchains
move_ratio[irot] = (2.0*ns_data.nchains) if ns_data.nbeads >= 2 else 0
move_ratio[idih] = 1.0*max(((ns_data.nbeads-3.0)*(ns_data.nchains),0))
move_ratio[ishear] = 3
move_ratio[istr] = 3

# print(move_ratio)

dof = 0


if move_ratio[0] != 0:
    dof+= 3*ns_data.nchains
if move_ratio[1] != 0:
    dof+= 2*ns_data.nchains



if parameters.previous_iterations == 0:

    print("Creating New Output File")

    energies_file.write(f"{ns_data.nwalkers} {1} {dof} {False} {ns_data.nchains} \n")
    energies_file.write(f"{ns_data.walklength} {ns_data.nbeads} {ns_data.nchains} {ns_data.nwalkers} \n")

energies_file.close()

###################################################################
active_box = ns_data.nwalkers+1


initialise_sim_cells(ns_data)


#random seeds
#np.random.seed(1)
#alk.random_set_random_seed(1)

###################################################################

#creating initial boxes/loading boxes from restart folder

parameters.configure_system()
    
    
#setting step sizes
alk.alkane_set_dr_max(0.65)
alk.alkane_set_dt_max(0.43)
alk.alkane_set_dh_max(0.4)
alk.alkane_set_dv_max(0.5)

ns_data.set_dshear_max(0.5)
ns_data.set_dstretch_max(0.5)



#moves_prob = np.cumsum(moves_ratio)/np.sum(moves_ratio)

dof_sum = np.sum(move_ratio)-7


#constructing dictionaries which contain initial volumes

if parameters.previous_iterations == 0:
    ns_data.perturb_initial_configs(move_ratio)
else:
    ns_data.load_volumes()
    
# main driver code
##############################################################

# mc_adjust_interval = ns_data.nwalkers//2 #for adjusting step sizes

snapshots = 100
vis_interval = max(1,parameters.total_iterations//snapshots)

# restart_interval = int(5e3)
# print_interval = int(1e2)

ns_data.set_intervals(vis_interval = vis_interval)

# for i in range(prev_lines, ns_iterations+prev_lines):
#     perform_ns_iter(ns_data, i, move_ratio, thread = 0)

print(parameters.previous_iterations,parameters.total_iterations)
perform_ns_run(ns_data,parameters.total_iterations, prev_iters=parameters.previous_iterations, move_ratio = move_ratio, processes=1, verbose = True)
        
      

energies_file.close()

write_configs_to_hdf(ns_data,parameters.previous_iterations+parameters.total_iterations)

#overlap check
ns_data.check_overlaps()

time_taken = ns_data.time_elapsed()

print(f"---------Time elapsed is {int(time_taken//60):0>2}:{time_taken%60:0>6.3f}---------")