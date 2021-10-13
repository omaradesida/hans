
from timeit import default_timer as timer

t0 = timer()

from NS_hsa import *

#parsing inputs

#parent_dir = os.getcwd()
parent_dir = "." #temporary until I rebuild hs_alkane to take longer filenames

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iterations", default = 2e5, type=float, help = "Number of iterations to run \n")
parser.add_argument("-t", "--time", type=float, help = "How much time has been allocated for the program to run")
parser.add_argument("-w", "--nwalkers", default = 500, type=int, help = "Number of walkers to use \n")
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




    
print("importing stuff")

#parsing inputs

parent_dir = os.getcwd()
parent_dir = "." #temporary until I rebuild hs_alkane to take longer filenames

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iterations", default = 2e5, type=float, help = "Number of iterations to run \n")
parser.add_argument("-w", "--nwalkers", default = 500, type=int, help = "Number of walkers to use \n")
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
    os.mkdir(f"{pwd}restarts/")
    
    

if from_restart:
    print("loading settings from prev run")
    restart_folder = args.restart_folder
    
    pwd = f"{parent_dir}/{restart_folder}"
    print(pwd)
    with open(f'{pwd}energies.txt') as f:
        dataline = f.readline().rstrip().split()
        walk_length = int(dataline[0])
        nbeads = int(dataline[1])
        nchains = int(dataline[2])
        nwalkers = int(dataline[3])

    
energies_filename = f"{pwd}energies.txt"
traj_filename = f"{pwd}traj.extxyz"



###################################################################
active_box = nwalkers+1

mdl.box_set_num_boxes(nwalkers+1) #nwalkers+2 if debugging
mdl.box_initialise()
mdl.box_set_pbc(1)
mdl.alkane_set_nchains(nchains) 
mdl.alkane_set_nbeads(nbeads)    
mdl.alkane_initialise()           
mdl.box_set_isotropic(1)
mdl.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
mdl.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate these

#random seeds
#np.random.seed(1)
#mdl.random_set_random_seed(1)

###################################################################

#creating initial boxes/loading boxes from restart folder

if not from_restart:
    max_vol_per_atom = 15
    #set cell size
    cell_matrix = 0.999*np.eye(3)*np.cbrt(nbeads*nchains*max_vol_per_atom)#*np.random.uniform(0,1)
    for ibox in range(1,nwalkers+1):
        mdl.box_set_cell(ibox,cell_matrix)
    populate_boxes(nwalkers, nchains)
    
else:
    mdl.io_read_xmol(f"{pwd}restarts/chain.xmol")
    
    
#setting step sizes
mdl.alkane_set_dr_max(0.65)
mdl.alkane_set_dt_max(0.43)
mdl.alkane_set_dh_max(0.4)
mdl.alkane_set_dv_max(0.5)


#setting params for MC_run
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']
ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5
move_ratio = np.zeros(6)
move_ratio[ivol] = 1
move_ratio[itrans] = 3.0*nchains
move_ratio[irot] = (2.0*nchains) if nbeads >= 2 else 0
move_ratio[idih] = max(((nbeads-3.0)*(nchains),0))
move_ratio[ishear] = 3
move_ratio[istr] = 3
#moves_prob = np.cumsum(moves_ratio)/np.sum(moves_ratio)


#constructing dictionaries which contain initial volumes

if not from_restart:
    volumes = shake_initial_configs(nwalkers, move_ratio)
else:
    volumes = {}
    for ibox in range(1,nwalkers+1):
        volumes[ibox] = mdl.box_compute_volume(ibox)
    
# main driver code

mc_adjust_interval = nwalkers//2 #for adjusting step sizes
low_acc_rate = 0.2
high_acc_rate = 0.5 

snapshots = 1000
vis_interval = ns_iterations//snapshots



energies_file = open(energies_filename, "a+")

if not from_restart:
    energies_file.write(f"{walk_length} {nbeads} {nchains} {nwalkers} \n")
    energies_file.write(f"{nwalkers} {1} {5*nchains} {False} {nchains} \n")
    prev_lines = 0
else:
    prev_lines = sum(1 for line in open(energies_filename)) - 2

    
    
restart_interval = int(5e3)

#traj = ase.io.write("nestedsampling.extxyz", mode="a", fmt="extxyz")


# f = IntProgress(min=0, max=ns_iterations) 
# display(f) # display the bar



for i in range(ns_iterations):
    index_max = max(volumes, key=volumes.get)
    volume_limit = volumes[index_max]
    clone = np.random.randint(1,nwalkers+1)
    
    if i%vis_interval == 0:
        largest_config = mk_ase_config(index_max,nbeads,nchains)
        largest_config.wrap()
        asewrite(traj_filename,largest_config, append = True)
#         traj.write(largest_config)


    
    if i%mc_adjust_interval == 0:
        adjust_dv(clone,active_box,low_acc_rate,high_acc_rate, volume_limit)
        adjust_dr(clone,active_box,low_acc_rate,high_acc_rate)
        if nbeads >= 2:
            adjust_dt(clone,active_box,low_acc_rate,high_acc_rate)
        if nbeads >= 4:
            adjust_dh(clone,active_box,low_acc_rate,high_acc_rate)
        adjust_dshear(clone,active_box,low_acc_rate,high_acc_rate)
        adjust_dstretch(clone,active_box,low_acc_rate,high_acc_rate)

    
    clone_walker(clone,active_box) #copies the ibox from first argument onto the second one.
    
    
    new_volume,_ = MC_run(walk_length, move_ratio, active_box,volume_limit)
    
    
    #removed_volumes.append(volumes[index_max])
    energies_file.write(f"{i+prev_lines} {volumes[index_max]} {volumes[index_max]} \n")
    volumes[index_max] = new_volume #replacing the highest volume walker
    clone_walker(active_box, index_max)
    if i%100 == 0:
        print(i,volume_limit)



    if i%restart_interval == 0:
        for j in range(1,nwalkers+2):
            restart_file = f"{pwd}restarts/chain.xmol.{j:04}"
            config = mk_ase_config(j,nbeads,nchains, scaling = 1)
            config.wrap()
            asewrite(restart_file, config, format = "extxyz",
                     comment = celltoxmolstring(config))
        t1 = timer()
        t_elapsed = t1-t0
        if t1-t_elapsed < 600:
            sys.exit()


        
energies_file.close()

for j in range(1,nwalkers+2):
    restart_file = f"{pwd}restarts/chain.xmol.{j:04}"
    config.wrap()
    config = mk_ase_config(j,nbeads,nchains, scaling = 1)
    asewrite(restart_file, config, format = "extxyz",
             comment = celltoxmolstring(config))

#overlap check
overlap_check = np.zeros(nwalkers)
for ibox in range(1,nwalkers+1):
    overlap_check[mdl.alkane_check_chain_overlap(ibox)]
if np.any(overlap_check):
    print("Overlaps Detected")
else:
    print("No Overlaps Present")
