import h5py
import NS_hsa as NS
from mpi4py import MPI

def read_hans_file(filename: str = "input.txt"):
    """Parse input from a file"""
    data = {}
    f = open(filename, "r")
    inputs=f.readlines()
    for line in inputs:
        key,value=line.split("=")
        data[key.strip()] = value.strip()
    f.close()
    data["nchains"] = int(data["nchains"])
    data["nbeads"] = int(data["nbeads"])
    data["nwalkers"] = int(data["nwalkers"])
    data["walklength"] = int(data["walklength"])
    data["iterations"] = float(data["iterations"])
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

def write_to_restart(filename,args,rank,current_iter=0,comm=None):
    f=h5py.File(filename,'w',driver='mpio',comm=MPI.COMM_WORLD)
    for i in args:
        f.attrs.create(i,args[i])
    for iwalker in range(1,args["nwalkers"]+1):
        groupname = f"walker_{iwalker:04d}"
        f.create_group(groupname)
    for i in range(rank*args["nwalkers"], (rank+1)*args["nwalkers"]):
        tmpgroup = f[f"walker_{iwalker:04d}"]

    
    return

