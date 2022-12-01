import h5py
from NesSa import MCNS as NS
from mpi4py import MPI
import argparse
import sys
import ase.io

def parse_args():

    """Parse command line arguments. Returns an empty Namespace containing the attributes."""

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--input_file",type=str,default = "input.txt", help =
                        "Which input file to use for the nested sampling simulation.\n")


    return parser.parse_args()

def read_hans_file():
    """Parse input from a file"""
    data = {}
    f = sys.stdin.read()
    if f == "":
        print("Error: No input file specified")
        sys.exit(1)
    inputs=f.split("\n")
    for line in inputs:
        if (not line.startswith('#') and line != ''):
            key,value=line.split("=")
            data[key.strip()] = value.strip()
    float_keys = ["bondlength","bondangle","min_angle","min_aspect_ratio", 
                "pressure", "upper_bound", "lower_bound", "time"]
    int_keys = ["nchains","nbeads","nwalkers","walklength","initial_walk"
                ,"analyse", "equil_iter", "main_iter", "index", "model_type"]
    for key in float_keys:
        if key in data:
            data[key] = float(data[key])
    for key in int_keys:
        if key in data:
            data[key] = int(data[key])
    if not "initial_walk" in data:
        data["initial_walk"] = 1000
    if not "analyse" in data:
        data["analyse"] = 0
    if not "bondlength" in data:
        data["bondlength"] = 0.4
    if not "bondangle" in data:
        data["bondangle"] = 109.7
    if not "min_angle" in data:
        data["min_angle"] = 45.0
    if  not "min_aspect_ratio" in data:
        data["min_aspect_ratio"] = 0.8
    if not "min_aspect_ratio" in data:
        data["min_aspect_ratio"] = 0.8
    if not "restart_file" in data:
        data["restart_file"] = "restart.hdf5"
    if not "model_type" in data:
        data["model_type"] = 4



    
    return data

def write_to_restart(args,comm,filename = "restart.hdf5",i=0):
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank==0:
        f = h5py.File(filename, "w")
        for j in args:
            f.attrs.create(j,args[j])
        f.attrs.create("prev_iters",i+1)
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
            config_list = comm.recv(source=j,tag = j)
            for iwalker in range(1,args["nwalkers"]+1):
                groupname = f"walker_{j}_{iwalker:04d}"
                tempgrp = f.create_group(groupname)
                coords = tempgrp.create_dataset("coordinates",(args["nbeads"]*args["nchains"],3),dtype="float64")
                unitcell = tempgrp.create_dataset("unitcell",(3,3),dtype="float64")

                unitcell[:] = config_list[iwalker-1].cell
                for ichain in range(args["nchains"]):
                    coords[:] = config_list[iwalker-1].positions
        f.close()
    else:
        config_list=[NS.mk_ase_config(ibox+1,args["nbeads"],args["nchains"],1.0) for ibox in range(args["nwalkers"])]
        comm.ssend(config_list,0,tag=rank)

def write_to_extxyz(args,ibox=1,filename="traj.extxyz", parallel = False):
    """Writes a single simulation box to file.
        Arguments:
            ibox: Simulation box to write. If none, the largest simulation box is used."""
    nbeads = args["nbeads"]
    nchains = args["nchains"]


    max_vol_config = NS.mk_ase_config(ibox,nbeads,nchains, scaling = 1.0)
    max_vol_config.wrap()

    NS.io.write(filename, max_vol_config, append = True, parallel=parallel)
    return

def restart_cleanup(args,traj_interval=100):
    old_vol_file = open("volumes.txt", "r+")
    lines  = old_vol_file.readlines()
    old_vol_file.close()
    if (len(lines) - 1) > args["prev_iters"]:
        print("Warning, restarting from an older file, data may be overwritten/deleted")
        sys.stdout.flush()
        lines = lines[:(args["prev_iters"]+1)]
        new_vol_file =  open("volumes.txt","w+")
        
        for line in lines:
            new_vol_file.write(line)
        sys.stdout.flush()
        new_vol_file.close()
        n_ase_images  = max(args["prev_iters"]//traj_interval,1)
        old_images = ase.io.read("traj.extxyz",f":{n_ase_images}", parallel = False)
        ase.io.write("traj.extxyz", old_images, parallel = False, append = False )
