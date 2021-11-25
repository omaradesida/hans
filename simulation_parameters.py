import os
import h5py
import numpy as np
import hs_alkane.alkane as alk
import cli

class SimulationParameters:

    def __init__(self, nwalkers, nchains, nbeads, walklength, directory, time, previous_iterations, total_iterations, processes, bondlength):
        self.nwalkers = nwalkers
        self.nchains = nchains
        self.nbeads = nbeads
        self.walklength = walklength
        self.directory = directory
        self.allotted_time = time
        self.previous_iterations = previous_iterations
        self.total_iterations = total_iterations
        self.processes = processes
        self.bondlength = bondlength

    @classmethod
    def parse_args(cls, parent_dir="."):
        args = cli.parse_args()
        for arg in vars(args):
            print (f"{arg:<15} {getattr(args, arg)}")


        if args.restart:
            parameters = SimulationParameters.from_restart(args, parent_dir)
        else:
            parameters = SimulationParameters.from_args(args, parent_dir)
        return parameters

    @classmethod
    def from_args(cls, args, parent_dir="."):
        nwalkers = args.nwalkers
        nchains = args.nchains
        nbeads = args.nbeads
        walklength = int(args.walklength)
        total_iterations = int(args.iterations)
        processes = int(args.processes)
        bondlength=args.bondlength
        
        dir_prefix = f"{parent_dir}/NS_{nchains}_{nbeads}mer.{nwalkers}.{walklength}"
        i_n = 1
        
        while os.path.exists(f"{dir_prefix}.{i_n}/"):
            i_n += 1
        
        directory = f"{dir_prefix}.{i_n}/"

        print(f"Creating directory {directory} to output results")
        os.mkdir(f"{directory}")

        return cls(nwalkers, nchains, nbeads, int(walklength), directory, args.time, 0, total_iterations, processes, bondlength)

    @classmethod
    def from_restart(cls, args, parent_dir="."):
        print("loading settings from prev run")
        restart_folder = args.restart_folder
        processes = int(args.processes)
        
        directory = f"{parent_dir}/{restart_folder}"

        if directory[-1] != "/":
            directory+="/"


        print(f"Attempting to restart from {directory}")

        f = h5py.File(f"{directory}restart.hdf5", "r")

        nbeads = int(f.attrs["nbeads"])
        nchains = int(f.attrs["nchains"])
        nwalkers = int(f.attrs["nwalkers"])
        previous_iterations = int(f.attrs["prev_iters"])
        walklength = int(f.attrs["sweeps"])
        total_iterations = int(args.iterations)
        bondlength = int(f.attrs['bondlength'])

        f.close()

        return cls(nwalkers, nchains, nbeads, walklength, directory, args.time, previous_iterations, total_iterations, processes, bondlength)

    def configure_system(self, max_vol_per_atom=15):
        if self.previous_iterations == 0:
            self.create_initial_configs(max_vol_per_atom)
        else:
            self.set_configs_from_hdf()

    def create_initial_configs(self, max_vol_per_atom = 15):
        cell_matrix = 0.999*np.eye(3)*np.cbrt(self.nbeads*self.nchains*max_vol_per_atom)#*np.random.uniform(0,1)
        for ibox in range(1,self.nwalkers+1):
            alk.box_set_cell(int(ibox),cell_matrix)
        self.populate_boxes()

    def populate_boxes(self):
        ncopy = self.nchains
        for ibox in range(1,self.nwalkers+1):
            for ichain in range(1,ncopy+1):
                rb_factor = 0
                alk.alkane_set_nchains(ichain)
                overlap_flag = 1
                while rb_factor == 0:
                    rb_factor, ifail = alk.alkane_grow_chain(ichain,int(ibox),1) 
                    if ifail != 0:
                        rb_factor = 0

    def set_configs_from_hdf(self):
        filename = f"{self.directory}restart.hdf5"
        f = h5py.File(filename,"r")

        nbeads = int(f.attrs["nbeads"])
        nchains = int(f.attrs["nchains"])
        nwalkers = int(f.attrs["nwalkers"])

        for iwalker in range(1,nwalkers+1):
            groupname = f"walker_{iwalker:04d}"
            cell = f[groupname]["unitcell"][:]
            alk.box_set_cell(iwalker,cell)
            new_coords = f[groupname]["coordinates"][:]
            for ichain in range(0,nchains):
                coords = alk.alkane_get_chain(ichain+1,iwalker)
                for ibead in range(nbeads):
                    coords[ibead] = new_coords[ichain*nbeads+ibead]

        f.close()


