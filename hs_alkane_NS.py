
from timeit import default_timer as timer

t0 = timer()

from NS_hsa import *

#parsing inputs

#parent_dir = os.getcwd()
parent_dir = "." #temporary until hs_alkane can take longer filenames

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iterations", default = 2e5, type=float, help = "Number of iterations to run \n")
parser.add_argument("-t", "--time", type=float, help = "How much time has been allocated for the program to run")
parser.add_argument("-w", "--nwalkers", default = 10, type=int, help = "Number of walkers to use \n")
parser.add_argument("-c", "--nchains", default = 32, type=int, help = "Number of chains in the simulation box \n")
parser.add_argument("-b", "--nbeads", default = 1, type=int, help = "Number of beads in each chain \n")
parser.add_argument("-l", "--walklength", default = 80, type=int, help = "Number of sweep for each random walk step \n")
parser.add_argument("-R", "--restart", action = "store_true", help = "Whether or not to restart from a previous attempt\n")
parser.add_argument("-f","--restart_folder",type=str, help =
                    "Configurations file to restart from. Energies and Trajectory file should have the same file naming \
                    convention. \n i.e if restart file is 'foo.restart' energies and trajectory file \
                    should be foo.energies and foo.extxyz")

args = parser.parse_args()
from_restart = args.restart
ns_iterations = int(args.iterations)
time = args.time



if not from_restart:
    nwalkers = args.nwalkers
    nchains = args.nchains
    nbeads = args.nbeads
    
    walk_length = int(args.walklength)
    
    dir_prefix = f"{parent_dir}/NS_{nchains}_{nbeads}mer.{nwalkers}.{walk_length}"
    i_n = 1
    
    while os.path.exists(f"{dir_prefix}.{i_n}/"):
        i_n += 1
            
#     traj_filename = f"{dir_prefix}.{i_n}.extxyz"
#     energies_filename = f"{dir_prefix}.{i_n}.energies"  
    
    pwd  = f"{dir_prefix}.{i_n}/"    
    os.mkdir(f"{pwd}")


if from_restart:
    print("loading settings from prev run")
    restart_folder = args.restart_folder
    
    pwd = f"{parent_dir}/{restart_folder}"
    print(pwd)
    with open(f'{pwd}energies.txt') as f:
        next(f)
        dataline = f.readline().rstrip().split()
        if len(dataline) != 4:
            raise ValueError("Wrong number of inputs in restart file")
        walk_length = int(dataline[0])
        nbeads = int(dataline[1])
        nchains = int(dataline[2])
        nwalkers = int(dataline[3])



ns_run = ns_info(nwalkers,nchains,nbeads,walk_length, time)

ns_run.set_directory(pwd)


###################################################################
active_box = ns_run.nwalkers+1

alk.box_set_num_boxes(ns_run.nwalkers+1) #nwalkers+2 if debugging
alk.box_initialise()
alk.box_set_pbc(1)
alk.alkane_set_nchains(ns_run.nchains) 
alk.alkane_set_nbeads(ns_run.nbeads)    
alk.alkane_initialise()           
alk.box_set_isotropic(1)
alk.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
alk.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate these

# ns_run.initialise_hsa()

#random seeds
#np.random.seed(1)
#alk.random_set_random_seed(1)

###################################################################

#creating initial boxes/loading boxes from restart folder

if not from_restart:

    create_initial_configs(ns_run)
    
else:
    alk.io_read_xmol(f"{pwd}restarts/chain.xmol")
    
    
#setting step sizes
alk.alkane_set_dr_max(0.65)
alk.alkane_set_dt_max(0.43)
alk.alkane_set_dh_max(0.4)
alk.alkane_set_dv_max(0.5)

ns_run.set_dshear_max(0.5)
ns_run.set_dstretch_max(0.5)


#setting params for MC_run
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']
ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5
move_ratio = np.zeros(6)
move_ratio[ivol] = 1
move_ratio[itrans] = 3.0*ns_run.nchains
move_ratio[irot] = (2.0*ns_run.nchains) if ns_run.nbeads >= 2 else 0
move_ratio[idih] = max(((ns_run.nbeads-3.0)*(ns_run.nchains),0))
move_ratio[ishear] = 3
move_ratio[istr] = 3
#moves_prob = np.cumsum(moves_ratio)/np.sum(moves_ratio)


#constructing dictionaries which contain initial volumes

if not from_restart:
    ns_run.perturb_initial_configs(move_ratio)
else:
    ns_run.load_volumes()
    
# main driver code

# mc_adjust_interval = ns_run.nwalkers//2 #for adjusting step sizes


# snapshots = 1000
# vis_interval = ns_iterations//snapshots

# restart_interval = int(5e3)
# print_interval = int(1e2)

ns_run.set_intervals()


#energies_file = open(ns_run.energies_filename, "a+")

if not from_restart:

    ns_run.energies_file.write(f"{ns_run.nwalkers} {1} {5*ns_run.nchains} {False} {ns_run.nchains} \n")
    ns_run.energies_file.write(f"{ns_run.sweeps_per_walk} {ns_run.nbeads} {ns_run.nchains} {ns_run.nwalkers} \n")
    prev_lines = 0
else:
    prev_lines = sum(1 for line in open(ns_run.energies_filename)) - 2

    

#traj = ase.io.write("nestedsampling.extxyz", mode="a", fmt="extxyz")


# f = IntProgress(min=0, max=ns_iterations) 
# display(f) # display the bar



for i in range(prev_lines, ns_iterations+prev_lines):
    perform_ns_iter(ns_run, i, move_ratio)

# for i in range(ns_iterations):
#     index_max = ns_run.max_vol_index()
#     volume_limit = ns_run.volumes[index_max]
#     clone = np.random.randint(1,nwalkers+1)
    
#     if i%vis_interval == 0:
#         ns_run.write_to_traj()
# #         traj.write(largest_config)
    
#     if i%mc_adjust_interval == 0:
#         adjust_mc_steps(ns_run, clone, active_box, volume_limit)

    
#     clone_walker(clone,active_box) #copies the ibox from first argument onto the second one.
    
    
#     new_volume,_ = MC_run(ns_run,walk_length, move_ratio, active_box,volume_limit)
    
    
#     #removed_volumes.append(volumes[index_max])
#     ns_run.energies_file.write(f"{i+prev_lines} {ns_run.volumes[index_max]} {ns_run.volumes[index_max]} \n")
#     ns_run.volumes[index_max] = new_volume #replacing the highest volume walker
#     clone_walker(active_box, index_max)
#     if i%print_interval == 0:
#         print(i,volume_limit)


#     if i%restart_interval == 0:
#         write_configs_to_hdf(ns_run,f"{pwd}restart.hdf5")
#         t1 = timer()
#         t_elapsed = t1-t0
#         if time-t_elapsed < 600:

#             print("Out of allocated time, writing to file and exiting")
#             sys.exit()


        
ns_run.energies_file.close()

write_configs_to_hdf(ns_run,f"{pwd}finalconfigs.hdf5")

#overlap check
ns_run.check_overlaps()