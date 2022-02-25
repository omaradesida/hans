import h5py
import NS_hsa as NS
from mpi4py import MPI

def read_hans_file(filename: str = "input.txt"):
    """Parse input from a file"""
    data = {}
    f = open(filename, "r")
    inputs=f.readlines()
    inputs=[line for line in inputs if line]
    for line in inputs:
        if not line.startswith('#'):
            print(line)
            key,value=line.split("=")
            data[key.strip()] = value.strip()
    f.close()
    data["nchains"] = int(data["nchains"])
    data["nbeads"] = int(data["nbeads"])
    data["nwalkers"] = int(data["nwalkers"])
    data["walklength"] = int(data["walklength"])
    data["iterations"] = float(data["iterations"])
    data["time"] = float(data["time"])
    if "bondlength" in data:
        data["bondlength"] = float(data["bondlength"])
    else:
        data["bondlength"] = 0.4
    if "bondangle" in data:
        data["bondangle"] = float(data["bondangle"])
    else:
        data["bondangle"] = 109.7
    if "min_angle" in data:
        data["min_angle"] = float(data["bondlength"])
    else:
        data["min_angle"] = 45.0
    if "min_aspect_ratio" in data:
        data["min_aspect_ratio"] = float(data["min_aspect_ratio"])
    else:
        data["min_aspect_ratio"] = 0.8
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


    max_vol_config = NS.mk_ase_config(ibox,nbeads,nchains)
    max_vol_config.wrap()

    NS.io.write(filename, max_vol_config, append = True, parallel=parallel)
    return