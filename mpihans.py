from timeit import default_timer as timer
import sys
from mpi4py import MPI
import NS_hsa as NS
from hans_io import read_hans_file
from numpy.random import MT19937
from numpy.random import RandomState, SeedSequence
import os
import ase
import h5py
import numpy as np
import cProfile

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def main():

    


    #constants for MC adjust stuff
    dv_max = 50.0 #max vol move allowed
    dr_max = 50.0 #max trans move allowed

    dv_min = 1e-10 #smallest moves allowed
    dr_min = 1e-10 
    dt_min = 1e-10
    dh_min = 1e-10
    dshr_min = 1e-5
    dstr_min = 1e-5

    lower_bound = 0.2 #mc adjust upper and lower bounds
    upper_bound = 0.5

    mc_adjust_wl = max(10//size,1)

    dstretch = 1.0  # initial stretch and shear move sizes
    dshear = 1.0

    #
    directory = None
    args = None
    if rank == 0:
        args = read_hans_file() #parsing input file
        args["nchains"] = int(args["nchains"])
        args["nbeads"] = int(args["nbeads"])
        args["nwalkers"] = int(args["nwalkers"])
        args["walklength"] = int(args["walklength"])
        args["bondlength"] = float(args["bondlength"])
        args["bondangle"] = float(args["bondangle"])
        args["iterations"] = float(args["iterations"])
        args["min_angle"] = float(args["min_angle"])
        args["min_aspect_ratio"] = float(args["min_aspect_ratio"])
        vol_file = open("volumes.txt","w+")

        if "directory" in args:
            dir_prefix = args["directory"]
        else:
            dir_prefix = f"NesSam_{args['nchains']}_{args['nbeads']}mer.{args['nwalkers']}.{args['walklength']}"
            i_n = 1

        while os.path.exists(f"{dir_prefix}.{i_n}/"):
            i_n += 1

        directory = f"{dir_prefix}.{i_n}/"

        for arg in args:
                print (f"{arg:<16} {args[arg]}")

        print(f"Creating directory {directory} to output results")
        os.mkdir(f"{directory}")

    args = comm.bcast(args,root=0) #broadcasting input arguments
    directory = comm.bcast(directory,root=0) #broadcasting input arguments

    move_ratio=NS.default_move_ratio(args) #generating a move ratio

    rs = RandomState(MT19937(SeedSequence(2451))) #random seed

    NS.initialise_sim_cells(args) #initialise data structures
    NS.create_initial_configs(args) #creating initial configs

    NS.alk.alkane_set_dr_max(0.65) #set step sizes
    NS.alk.alkane_set_dt_max(0.43)
    NS.alk.alkane_set_dh_max(0.4)
    NS.alk.alkane_set_dv_max(0.5)

    NS.perturb_initial_configs(args,move_ratio, 50) #generating random box sizes

    vols=[NS.alk.box_compute_volume(i) for i in range(1,args["nwalkers"]+1)]

    mc_adjust_interval = (args["nwalkers"]*size)//2 #ns_adjust interval steps, same as pymatnest
    total_rate = np.zeros(6)

    dof = 0
    if move_ratio[0] != 0:
        dof+= args["nchains"]
    if move_ratio[1] != 0:
        if args["nbeads"] <= 2:
            dof+= 2*args["nchains"]
        else:
            dof+= 3*args["nchains"] #kinetic degrees of freedom for ns_analyse

    f = None
    if rank == 0:
        f = open(f"{directory}volumes.txt","w+")
        f.write(f'{args["nwalkers"]} {1} {dof} {False} {args["nchains"]} \n')

    #comm.Barrier()
    t0 = timer()
    i = 0
    for i in range(int(args["iterations"])):
        local_max = max(vols)
        local_max_index = [rank,vols.index(local_max)]


        vol_max,vol_max_index = comm.allreduce([local_max,local_max_index],op=MPI.MAXLOC)


        walker_to_clone = None
        if rank == 0:
            walker_to_clone=divmod(np.random.randint(args["nwalkers"]*size),args["nwalkers"])
            f.write(f"{i} {vol_max:.13f} {vol_max:.13f} \n")

        walker_to_clone=comm.bcast(walker_to_clone,root=0)
        comm.Barrier()



        config_to_clone=None
        if rank==walker_to_clone[0]:
            #print(i,walker_to_clone[0],vol_max_index[0])        
            config_to_clone = NS.mk_ase_config(walker_to_clone[1]+1,args["nbeads"],args["nchains"],scaling=1.0)
        config_to_clone = comm.bcast(config_to_clone, root = walker_to_clone[0])#,dest=vol_max_index[0])


        if rank == vol_max_index[0]:
            if i%100 == 0:
                NS.write_to_extxyz(args,vol_max_index[1]+1, filename=f"{directory}traj.extxyz")
                print(i)
            active_walker = vol_max_index[1]
            NS.import_ase_to_ibox(config_to_clone,active_walker+1,args)
        else:
            active_walker = np.random.randint(args["nwalkers"])

        comm.Barrier()
        new_vol,_ = NS.MC_run_2(args,args["walklength"], move_ratio,active_walker+1, volume_limit=vol_max,
                                    min_ar=args["min_aspect_ratio"], min_ang= args["min_angle"],
                                    dshear = dshear, dstretch = dstretch)

        vols[active_walker] = new_vol


        if i%mc_adjust_interval == 0:
            #print(NS.alk.alkane_get_dv_max(),NS.alk.alkane_get_dr_max(),dshear,dstretch, rank)
            rate = np.zeros(6)
            avg_rate = np.zeros_like(rate)
            mc_box = np.random.randint(args["nwalkers"])

            for i in range(6):
                move_ratio_matrix = np.eye(6)
                backup = NS.mk_ase_config(mc_box+1,args["nbeads"],args["nchains"],scaling=1)
                if move_ratio[i] != 0:
                    rate += NS.MC_run_2(args,mc_adjust_wl, move_ratio_matrix[i],mc_box+1,vol_max,dshear=dshear, dstretch=dstretch)[1]
                    NS.import_ase_to_ibox(backup,mc_box+1,args)
            comm.Allreduce(rate,avg_rate,op=MPI.SUM)
            avg_rate = avg_rate/size
            if move_ratio[0] != 0:
                if avg_rate[0] < lower_bound:
                    NS.alk.alkane_set_dv_max(max(0.5*NS.alk.alkane_get_dv_max(),dv_min))
                elif avg_rate[0] > upper_bound:
                    NS.alk.alkane_set_dv_max(min(2.0*NS.alk.alkane_get_dv_max(),dv_max))
            if move_ratio[1] != 0:
                if avg_rate[1] < lower_bound:
                    NS.alk.alkane_set_dr_max(max(0.5*NS.alk.alkane_get_dr_max(),dr_min))
                elif avg_rate[1] > upper_bound:
                    NS.alk.alkane_set_dr_max(min(2.0*NS.alk.alkane_get_dr_max(),dr_max))
            if move_ratio[2] != 0:
                if avg_rate[2] < lower_bound:
                    NS.alk.alkane_set_dt_max(max(0.5*NS.alk.alkane_get_dt_max(),dt_min))
                elif avg_rate[2] > upper_bound:
                    NS.alk.alkane_set_dt_max(2.0*NS.alk.alkane_get_dt_max())
            if move_ratio[3] != 0:
                if avg_rate[3] < lower_bound:
                    NS.alk.alkane_set_dh_max(max(0.5*NS.alk.alkane_get_dh_max(),dh_min))
                elif avg_rate[3] > upper_bound:
                    NS.alk.alkane_set_dh_max(2.0*NS.alk.alkane_get_dh_max())
            if move_ratio[4] != 0:
                if avg_rate[4] < lower_bound:
                    dshear = max(0.5*dshear,dshr_min)
                elif avg_rate[4] > upper_bound:
                    dshear = 2.0*dshear
            if move_ratio[5] != 0:
                if avg_rate[5] < lower_bound:
                    dstretch = max(0.5*dstretch,dstr_min)
                elif avg_rate[5] > upper_bound:
                    dstretch  = 2.0*dstretch
            #print(NS.alk.alkane_get_dv_max(),NS.alk.alkane_get_dr_max(),dshear,dstretch, rank)
            #print(rate, rank)


        # #adjusting step_sizes
        # if i%mc_adjust_interval == 0:_,rate = NS.MC_run_2(args,mc_adjust_wl, [1,0,0,0,0,0],mc_box+1,vol_max)
        #     avgrate = np.zeros_like(total_rate)
        #     comm.Allreduce(total_rate,avgrate, op = MPI.SUM)
        #     avgrate = avgrate/(size*mc_adjust_interval)
        #     total_rate = np.zeros(6)
        # #comm.Allreduce(rate,avgrate,op=MPI.SUM)
        # #avgrate = avgrate/size
        # #adjust step sizes
        #     if move_ratio[0] != 0:
        #         if avgrate[0] < lower_bound:
        #             NS.alk.alkane_set_dv_max(max(0.5*NS.alk.alkane_get_dv_max(),dv_min))
        #         elif avgrate[0] > upper_bound:
        #             NS.alk.alkane_set_dv_max(min(2.0*NS.alk.alkane_get_dv_max(),dv_max))
        #     if move_ratio[1] != 0:
        #         if avgrate[1] < lower_bound:
        #             NS.alk.alkane_set_dr_max(max(0.5*NS.alk.alkane_get_dr_max(),dr_min))
        #         elif avgrate[1] > upper_bound:
        #             NS.alk.alkane_set_dr_max(min(2.0*NS.alk.alkane_get_dr_max(),dr_max))
        #     if move_ratio[2] != 0:
        #         if avgrate[2] < lower_bound:
        #             NS.alk.alkane_set_dt_max(max(0.5*NS.alk.alkane_get_dt_max(),dt_min))
        #         elif avgrate[2] > upper_bound:
        #             NS.alk.alkane_set_dt_max(2.0*NS.alk.alkane_get_dt_max())
        #     if move_ratio[3] != 0:
        #         if avgrate[3] < lower_bound:
        #             NS.alk.alkane_set_dh_max(min(0.5*NS.alk.alkane_get_dh_max(),dh_min))
        #         elif avgrate[i3] > upper_bound:
        #             NS.alk.alkane_set_dh_max(2.0*NS.alk.alkane_get_dh_max())
        #     if move_ratio[4] != 0:
        #         if avgrate[4] < lower_bound:
        #             dshear = 0.5*dshear
        #         elif avgrate[4] > upper_bound:
        #             dshear = max((2.0*dshear,dshr_min))
        #     if move_ratio[5] != 0:
        #         if avgrate[5] < lower_bound:
        #             dstretch = 0.5*dstretch
        #         elif avgrate[5] > upper_bound:
        #             dstretch  = max((2.0*dstretch,dstr_min))
            

    #comm.Barrier()

    if rank == 0:
        f.close()

    t2 = timer()
    if rank == 0:
        print("NS RUN TIME TAKEN =", t2-t0)            
        print("writing restart")
    if rank ==0:
        f = h5py.File("restart.hdf5", "w")
        for j in args:
            f.attrs.create(j,args[j])
        f.attrs.create("prev_iters",i)
        for iwalker in range(1,args["nwalkers"]+1):
            groupname = f"walker_{rank}_{iwalker:04d}"
            tempgrp = f.create_group(groupname)
            coords = tempgrp.create_dataset("coordinates",(args["nbeads"]*args["nchains"],3),dtype="float64")
            unitcell = tempgrp.create_dataset("unitcell",(3,3),dtype="float64")

            unitcell[:] = NS.alk.box_get_cell(iwalker)
            for ichain in range(args["nchains"]):
                chain = NS.alk.alkane_get_chain(ichain+1,iwalker)
                coords[ichain*args["nbeads"]:ichain*args["nbeads"]+(args["nbeads"]), :] = chain
        for j in range(1,size):
            print(j, " we start here")
            config_list = comm.recv(source=j,tag = j)
            for iwalker in range(1,args["nwalkers"]+1):
                groupname = f"walker_{j}_{iwalker:04d}"
                tempgrp = f.create_group(groupname)
                coords = tempgrp.create_dataset("coordinates",(args["nbeads"]*args["nchains"],3),dtype="float64")
            unitcell = tempgrp.create_dataset("unitcell",(3,3),dtype="float64")
            print(j, " We done here")
        f.close()
    else:
        config_list=[NS.mk_ase_config(ibox+1,args["nbeads"],args["nchains"],1.0) for ibox in range(args["nwalkers"])]
        comm.ssend(config_list,0,tag=rank)

    sys.stdout.flush()
    NS.alk.alkane_destroy()
    NS.alk.box_destroy()
    comm.Barrier()

    # print(NS.alk.alkane_get_dv_max(), "dv")
    # print(NS.alk.alkane_get_dr_max(), "dr")
    # print(NS.alk.alkane_get_dt_max(), "dt")
    # print(NS.alk.alkane_get_dh_max(), "dh")
if __name__ == "__main__":
    # prof = cProfile.Profile()
    # prof.enable()
    main()
    # prof.disable()
    # prof.dump_stats(f"mpihans{rank}.profile")
    print("All good G")

