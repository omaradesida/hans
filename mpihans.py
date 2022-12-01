from timeit import default_timer as timer
import sys
from mpi4py import MPI
from NesSa import MCNS as NS
from NesSa import NSio
#from numpy.random import MT19937
#from numpy.random import RandomState, SeedSequence
import os
import ase.io
import h5py
import numpy as np
from NesSa import ns_analyse_main

#import signal


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
t0 = MPI.Wtime()

def main(SimParams):
    #constants for MC adjust stuff
    dv_max = 50.0 #max vol move allowed
    dr_max = 50.0 #max trans move allowed

    dv_min = 1e-10 #smallest moves allowed
    dr_min = 1e-10 
    dt_min = 1e-10
    dh_min = 1e-10
    dshr_min = 1e-5
    dstr_min = 1e-5
    min_dstep = np.array([dv_min,dr_min,dt_min,dh_min,dshr_min,dstr_min])


    mc_aargsdjust_wl = max(10//size,1)

    traj_interval = 1000
    
    #############################################################################
    directory = None
    from_restart = None
    if rank == 0:
        #cl_args = hans_io.parse_args()
        if "from_restart" in SimParams and int(SimParams["from_restart"]):
                print(f"Attempting to restart from {SimParams['directory']}")
                from_restart = True
                directory = SimParams["directory"]
        else:
            from_restart = False
            if "directory" in SimParams:
                if SimParams["directory"] == "gen_prefix":
                    print("Generating folder(s) with appropriate prefix to for NS runs")
                    dir_prefix = f"NeSa_{SimParams['nchains']}_{SimParams['nbeads']}mer.{SimParams['nwalkers']*size}.{SimParams['walklength']*size}"
                    i_n = 1
                    while os.path.exists(f"{dir_prefix}.{i_n}/"):
                        i_n += 1
                    directory = f"{dir_prefix}.{i_n}/"
                else:
                    directory = SimParams["directory"]
                    print(f"Using directory {directory} specified in input file.")
                os.mkdir(f"{directory}")
            else:
                print("No directory specified, using default")
                directory = "./"


    directory = comm.bcast(directory,root=0)
    os.chdir(f"{directory}")

    SimParams = comm.bcast(SimParams,root=0)
    from_restart = comm.bcast(from_restart,root=0)

    if not "prev_iters" in SimParams:
        SimParams["prev_iters"] = 0

    if from_restart:
        f = h5py.File(SimParams["restart_file"], "r")

        for i in f.attrs:
            if i != "restart_file" and i != "time":
                SimParams[i] = f.attrs[i]
                if isinstance(SimParams[i],np.floating):
                    SimParams[i] = float(SimParams[i])
                elif isinstance(SimParams[i],np.integer):
                    SimParams[i] = int(SimParams[i])
                else:
                    continue
        f.close()
        if rank == 0:
                NSio.restart_cleanup(SimParams,traj_interval)
            # pass
            
        comm.Barrier()   
        
    else:
        if rank == 0:
            print("New Run")
            print("Files output to:",directory)
        
    if rank == 0:
        for arg in SimParams:
                print (f"{arg:<16} {SimParams[arg]}")


    move_ratio=NS.default_move_ratio(SimParams) #generating a move ratio
    SimParams["move_ratio"] = move_ratio #writing here so it gets written to restart

    # rs = RandomState(MT19937(SeedSequence(2451))) #random seed

    if rank:
        quiet = 1
    else:
        quiet = 0

    NS.initialise_sim_cells(SimParams,quiet = quiet) #initialise data structure

    if "dv_max" in SimParams:
        NS.alk.alkane_set_dv_max(float(SimParams["dv_max"])) #set step sizes
    else:
        NS.alk.alkane_set_dv_max(2.0)
    if "dr_max" in SimParams:
        NS.alk.alkane_set_dr_max(float(SimParams["dr_max"])) #set step sizes
    else:
        NS.alk.alkane_set_dr_max(2.0)
    if "dt_max" in SimParams:
        NS.alk.alkane_set_dt_max(float(SimParams["dt_max"])) #set step sizes
    else:
        NS.alk.alkane_set_dt_max(0.43)
    if "dh_max" in SimParams:
        NS.alk.alkane_set_dh_max(float(SimParams["dh_max"])) #set step sizes
    else:
        NS.alk.alkane_set_dh_max(0.4)


    if "dstretch" in SimParams:
        dstretch = SimParams["dstretch"] #set step sizes
    else:
        dstretch = 1.0

    if "dshear" in SimParams:
        dstretch = SimParams["dshear"] #set step sizes
    else:
        dshear = 1.0

    mc_adjust_wl = 10

    if not from_restart:
        if "initial_config" in SimParams:
            if rank == 0:
                print("Loading initial config")
            
            try:
                for i in range(SimParams["nwalkers"]):
                    initial_config = ase.io.read(f'../{SimParams["initial_config"]}')
                    NS.import_ase_to_ibox(initial_config,i+1,SimParams)
            except OSError:
                print(f"Cannot locate {SimParams['initial_config']} in {os.getcwd()}")
                sys.exit(1)
            if 'initial_config' in globals():
                assert(initial_config.get_number_of_atoms() == SimParams["nwalkers"]*SimParams["nbeads"]), "Initial config has wrong number of atoms"

        else:
            NS.create_initial_configs(SimParams) #creating initial configs
        NS.perturb_initial_configs(SimParams,move_ratio, SimParams["initial_walk"]) #random walk helps to distribute box sizes.
    else: # load from restart
        f = h5py.File(SimParams["restart_file"], "r")
        dshear = f.attrs["dshear"]
        dstretch = f.attrs["dstretch"]

        for iwalker in range(1,SimParams["nwalkers"]+1):
            groupname = f"walker_{rank}_{iwalker:04d}"
            cell = f[groupname]["unitcell"][:]
            NS.alk.box_set_cell(iwalker,cell)
            new_coords = f[groupname]["coordinates"][:]
            for ichain in range(0,SimParams["nchains"]):
                coords = NS.alk.alkane_get_chain(ichain+1,iwalker)
                for ibead in range(SimParams["nbeads"]):
                    coords[ibead] = new_coords[ichain*SimParams["nbeads"]+ibead]

        f.close()

    vols=[NS.alk.box_compute_volume(i) for i in range(1,SimParams["nwalkers"]+1)]

    mc_adjust_interval = max((SimParams["nwalkers"]*size)//2,1) #ns_adjust interval steps, same as pymatnest


#calculating degrees of freedom
    dof = 0
    if move_ratio[0] != 0:
        dof+= SimParams["nchains"]
    if move_ratio[1] != 0:
        if SimParams["nbeads"] <= 2:
            dof+= 2*SimParams["nchains"]
        else:
            dof+= 3*SimParams["nchains"] #kinetic degrees of freedom for ns_analyse

    f = None
    if rank == 0:
        f = open(f"volumes.txt","a+")
        if not from_restart:
            f.write(f'{SimParams["nwalkers"]*size} {1} {dof} {False} {SimParams["nchains"]} \n')
    sys.stdout.flush()
#######################################################################################
# NESTED SAMPLING LOOP                                                                #
#######################################################################################
    ns_t0 = timer()
    interrupted = False
    #signal.signal(signal.SIGTERM, NS.signal_handler)
    for i in range(SimParams["prev_iters"],SimParams["prev_iters"]+int(SimParams["iterations"])):
        local_max = max(vols)
        local_max_index = [rank,vols.index(local_max)]


        vol_max,vol_max_index = comm.allreduce([local_max,local_max_index],op=MPI.MAXLOC)


        walker_to_clone = None
        if rank == 0:
            walker_to_clone=divmod(np.random.randint(SimParams["nwalkers"]*size),SimParams["nwalkers"])
            f.write(f"{i} {vol_max:.13f} {vol_max:.13f} \n")

        walker_to_clone=comm.bcast(walker_to_clone,root=0)

        config_to_clone=None
        if rank==walker_to_clone[0]:
            #print(i,walker_to_clone[0],vol_max_index[0])        
            config_to_clone = NS.mk_ase_config(walker_to_clone[1]+1,SimParams["nbeads"],SimParams["nchains"],scaling=1.0)
        config_to_clone = comm.bcast(config_to_clone, root = walker_to_clone[0])#,dest=vol_max_index[0])


        if rank == vol_max_index[0]:
            if i%traj_interval == 0:
                NSio.write_to_extxyz(SimParams,vol_max_index[1]+1, filename=f"traj.extxyz")
                #print(i, vol_max)
            active_walker = vol_max_index[1]
            NS.import_ase_to_ibox(config_to_clone,active_walker+1,SimParams)
        else:
            active_walker = np.random.randint(SimParams["nwalkers"])

        new_vol,_ = NS.MC_run(SimParams,SimParams["walklength"], move_ratio,active_walker+1, volume_limit=vol_max,
                                    min_ar=SimParams["min_aspect_ratio"], min_ang= SimParams["min_angle"],
                                    dshear = dshear, dstretch = dstretch)

        vols[active_walker] = new_vol


        if i%mc_adjust_interval == 0:
            r, dshear,dstretch = NS.adjust_mc_steps(SimParams,comm,move_ratio,vol_max,walklength = mc_adjust_wl, 
                      min_dstep=min_dstep, dv_max=dv_max,dr_max=dr_max,dshear = dshear, dstretch = dstretch)
            #Adjusting length of step sizes based on trial acceptance rates.
            if rank == 0:
                print(i,vol_max,r)

        ####restart handler####

        t1 = MPI.Wtime()
        if rank == 0:            
            if (SimParams["time"] - (t1-t0)) < 300.0: 
                interrupted == True

        interrupted = comm.bcast(interrupted,root=0)
        if interrupted:
            if rank == 0:
                f.close()
                print("Out of allocated time, writing to file and exiting")
            break
        if (i+1) % 50000 ==0:
            if rank==0:
                if os.path.exists("restart_backup.hdf5"):
                    os.remove("restart_backup.hdf5")
                if os.path.exists("restart.hdf5"):
                    os.rename("restart.hdf5","restart_backup.hdf5")
            NSio.write_to_restart(SimParams,comm,filename = "restart.hdf5",i=i, dshear = dshear, dstretch = dstretch)
            sys.stdout.flush()
            if rank ==0:
                print("wrote to restart")
            # try:
            #     os.remove(f"restart.{i-100000}.hdf5")
            # except:
            #     pass

#######################################################################################
# END NESTED SAMPLING LOOP                                                            #
#######################################################################################
    if rank == 0:
        f.close()

    ns_t1 = timer()
    if rank == 0:
        print("NS RUN TIME TAKEN =", ns_t1-ns_t0)            
        print("writing restart")

        if os.path.exists("restart_backup.hdf5"):
            os.remove("restart_backup.hdf5")
        if os.path.exists("restart.hdf5"):
            os.rename("restart.hdf5","restart_backup.hdf5")
    NSio.write_to_restart(SimParams,comm,filename = "restart.hdf5",i=i,dshear=dshear, dstretch = dstretch)

    sys.stdout.flush()
    NS.alk.alkane_destroy()
    NS.alk.box_destroy()

    comm.Barrier()


###analysis####
    if SimParams["analyse"]:
        from contextlib import redirect_stdout
        print("running ns_analyse")
        nsa_args = {"T_min": 0.01,
                "dT": 0.001,
                "n_T": 1000,
                "kB": 1.0,
                "accurate_sum": False,
                "verbose":False,
                "profile":False,
                "skip":0,
                "line_end":None,
                "interval":1,                     
                "delta_pressure":0.0,
                "dump_terms":-1.0,
                "entropy":False,
                "quiet":False,
                "kolmogorov_smirnov_test":False,
                "files":"volumes.txt"}
        with open('ns_analyse.out', 'w') as f:
            with redirect_stdout(f):
                ns_analyse_main.analyse(nsa_args)
    

if __name__ == "__main__":
    SimParams = NSio.read_hans_file() #parsing input file
    if "profile" in SimParams:
        if SimParams["profile"]:
            import cProfile
            prof = cProfile.Profile()
            prof.enable()
    main(SimParams)
    if "profile" in SimParams:
        if SimParams["profile"]:
            prof.disable()
            prof.dump_stats(f"mpihans{rank}.profile")
    exit()

