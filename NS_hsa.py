# initial params/setting up empty boxes
import hs_alkane.alkane as alk
import numpy as np
import sys
import copy
from ase import io
from timeit import default_timer as timer
import multiprocessing as mp
from functools import partial
from itertools import repeat


#for reading/writing
import h5py



#for converting hs_alkane boxes to ASE atoms objects
from ase import Atoms
from ase.visualize import view
# from ase.io import write as asewrite



#Setting constants 

#Take this out at some point, its bad
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']

ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5

shear_step_max = 0.5
stretch_step_max = 0.5



class ns_info:

    def __init__(self, parameters, from_restart = False):
        self.parameters = parameters

        self.shear_step_max = float(0.5)
        self.stretch_step_max = float(0.5)

        self.low_acc_rate = float(0.2)
        self.high_acc_rate = float(0.5)
        
        self.iter_ = int(0)

        self.from_restart = from_restart

        self.start_time = timer()

        self.active_box = parameters.nwalkers+1

    # def __getattr__(self, item):
    #     return getattr(self.parameters, item)
    
    def time_elapsed(self):
        current_time = timer()
        self.time_elapsed_ = current_time - self.start_time
        return self.time_elapsed_

    def time_remaining(self):
        time_remaining_ = self.parameters.allotted_time - self.time_elapsed()

        return time_remaining_

    




    def set_dshear_max(self,dshear_max):
        self.shear_step_max = dshear_max

    def set_dstretch_max(self,dstretch_max):
        self.stretch_step_max = dstretch_max

    def set_acc_rate_range(self,acc_range):
        if len(acc_range) != 2:
            raise IndexError("'acc_range' should have two values in it")
        
        self.low_acc_rate = acc_range[0]
        self.high_acc_rate = acc_range[1]

    def set_intervals(self,mc_adjust_interval = None,vis_interval = None, 
                        restart_interval = int(5e3), print_interval = 100):
        if mc_adjust_interval is None:
            mc_adjust_interval = max(1,self.parameters.nwalkers//2)

        self.mc_adjust_interval = mc_adjust_interval
        
        if vis_interval is None:
            vis_interval = mc_adjust_interval

        self.vis_interval = vis_interval

        self.restart_interval = restart_interval

        self.print_interval = print_interval
        



    def check_overlaps(self):
        overlap_check = np.zeros(self.parameters.nwalkers)
        for ibox in range(1,self.parameters.nwalkers+1):
            overlap_check[alk.alkane_check_chain_overlap(int(ibox))]
        if np.any(overlap_check):
            print("Overlaps Detected")
            return 1
        else:
            print("No Overlaps Present")
            return 0

    def perturb_initial_configs(self, move_ratio, walk_length = 20):
    
        """ Runs a number of Monte Carlo steps on every simulation box, using the move_ratio assigned to it,
        Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume"""

        nwalkers = self.parameters.nwalkers

        self.volumes = {}
        #start_volumes = []
        for ibox in range(1,nwalkers+1):
            self.volumes[ibox], _ = MC_run(self, walk_length, move_ratio, ibox)

        #overlap check
        self.check_overlaps()

        return self.volumes


    def set_directory(self,path="./"):

        self.pwd = path

        self.restart_filename = f"{path}restart.hdf5"
        self.traj_filename = f"{path}traj.extxyz"
        self.energies_filename = f"{path}energies.txt"

        # self.energies_file = open(self.energies_filename, "a+")

        


            


    def load_volumes(self):
        self.volumes = {}
        for ibox in range(1,self.parameters.nwalkers+1):
            self.volumes[ibox] = alk.box_compute_volume(int(ibox))

        #overlap check
        self.check_overlaps()

        return self.volumes

    def write_to_extxyz(self, ibox=None):

        if ibox is None:
            ibox = self.max_vol_index()

        
        max_vol_config = mk_ase_config(ibox,self.parameters.nbeads,self.parameters.nchains)
        max_vol_config.wrap()

        io.write(self.traj_filename, max_vol_config, append = True)

    def max_vol_index(self):
        """ Returns index of largest simulation box"""
        self.max_vol_index_ = max(self.volumes, key=self.volumes.get)
        return self.max_vol_index_

    def write_all_to_extxyz(self, filename = "dump.extxyz"):
        for i in range(1,self.parameters.nwalkers+1):
            atoms = mk_ase_config(i,self.parameters.nbeads,self.parameters.nchains)
            atoms.wrap()
            io.write(filename, atoms, append = True)


def mk_ase_config(ibox, Nbeads, Nchains, scaling = 3.75):
    '''Uses the current state of the alkane model to construct an ASE atoms object'''
    
    # Create and populate ASE object
    model_positions = np.empty([Nchains*Nbeads, 3])
    cell_vectors = alk.box_get_cell(int(ibox))
   
    for ichain in range(0, Nchains):
        model_positions[Nbeads*ichain:Nbeads*ichain+Nbeads] = alk.alkane_get_chain(ichain+1, int(ibox))
    
    confstring = "C"+str(Nbeads*Nchains)
    
    box_config = Atoms(confstring, positions=model_positions*scaling, pbc=True, cell=cell_vectors*scaling)


    return box_config  # Returns ASE atom object

def vis_chains(vis_config, nbeads, nchains):
    '''Takes an ASE atoms object or list thereof and creates a customised ngl viewer 
       with appropriate settings for our bulk alkane chain models'''

    met = 0.35
    rad = 1.0
    
    colours = ['#DDDDDD']#['#FF1111','#FFAAAA', '#DDDDDD', '#1111FF', '#AAAAFF']
    ncols = len(colours)
    
    try:
            
        vis_config.wrap()
    except:
        pass
    
    sel=list()
    for icol in range(ncols):
        sel.append(list())
    
    # Create lists for each colour
    for ichain in range(0, nchains):
        icol = ichain%ncols
        for ibead in range(nbeads):
            iatom = ichain*nbeads + ibead
            sel[icol].append(iatom)
            
    v = view(vis_config, viewer='ngl')                   
    v.view.clear_representations()
    v.view.add_representation('unitcell', color='#000000')
    
    for icol in range(ncols):
        v.view.add_representation('ball+stick', selection=sel[icol], color=colours[icol], radius=rad, metalness=met)

    return v


def min_aspect_ratio(ibox):
    vol = alk.box_compute_volume(int(ibox))
    cell = alk.box_get_cell(int(ibox)).copy()
    test_dist = []
    #returns the shortest distance between two parallel faces, scaled such that the cell has a volume of 1
    min_aspect_ratio = sys.float_info.max
    
    for i in range(3):
        vi = cell[i,:]
        vnorm_hat = np.cross(cell[(i+1)%3,:],cell[(i+2)%3,:])
        vnorm_hat = vnorm_hat/(np.sqrt(np.dot(vnorm_hat,vnorm_hat)))
        min_aspect_ratio = min(min_aspect_ratio, abs(np.dot(vnorm_hat,vi)))
        test_dist.append(abs(np.dot(vnorm_hat,vi)))
        
    return min_aspect_ratio/np.cbrt(vol)


def min_angle(ibox):
    cell = alk.box_get_cell(int(ibox)).copy()
    min_angle = sys.float_info.max
    
    for i in range(3):
        vc1 = cell[(i+1)%3,:]
        vc2 = cell[(i+2)%3,:]
        dot_prod = np.dot(vc1,vc2)/np.sqrt((np.dot(vc1,vc1))*(np.dot(vc2,vc2)))
        if dot_prod < 0:
            dot_prod *= -1
        vec_angle = np.arccos(dot_prod)
        min_angle = min(min_angle,vec_angle)
        
    return min_angle


def box_shear_step(ibox, ns_data, aspect_ratio_limit = 0.8, angle_limit = 55):

    step_size = ns_data.shear_step_max
    # pick random vector for shear direction
    #np.random.seed(10)
    rnd_vec_ind = np.random.randint(0, 3)
    # turn other two into orthonormal pair
    #should I have a pair or should I have a single vector.
    other_vec_ind = list(range(3))
    other_vec_ind.remove(rnd_vec_ind)
    orig_cell = copy.deepcopy(alk.box_get_cell(int(ibox)))
    orig_cell_copy = copy.deepcopy(orig_cell)

    
    v1 = orig_cell_copy[other_vec_ind[0],:]
    v2 = orig_cell_copy[other_vec_ind[1],:]  

    
    v1 /= np.sqrt(np.dot(v1,v1))
    v2 -= v1*np.dot(v1,v2) 
    v2 /= np.sqrt(np.dot(v2,v2))

    if np.isnan(np.sum(v1)) or np.isnan(np.sum(v2)):
        print(v1,v2)
        print(orig_cell)
        sys.exit()



    
    # pick random magnitudes
    rv1 = np.random.uniform(-step_size, step_size)
    rv2 = np.random.uniform(-step_size, step_size)
#     rv1 = 2
#     rv2 = 0

    # create new cell and transformation matrix (matrix is additive)
    new_cell = orig_cell.copy()

    new_cell[rnd_vec_ind,:] += rv1*v1 + rv2*v2

    delta_H = new_cell - orig_cell
    
    #reject due to poor aspect ratio

    

    
    alk.alkane_change_box(int(ibox),delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    
    angle_limit_rad = angle_limit*np.pi/180
    
    if min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif alk.alkane_check_chain_overlap(int(ibox)):
        boltz = 0   
    else:
        boltz = 1
    
    
    #bake rejection  due to shape in here as it becomes easier to fit with the rest of the code
    
    
    return boltz, delta_H

def box_stretch_step(ibox,ns_data, aspect_ratio_limit = 0.8, angle_limit = 55):
    step_size = ns_data.stretch_step_max
    cell = alk.box_get_cell(int(ibox))
    new_cell = cell.copy()
    rnd_v1_ind = np.random.randint(0, 3)
    rnd_v2_ind = np.random.randint(0, 3)
    if rnd_v1_ind == rnd_v2_ind:
        rnd_v2_ind = (rnd_v2_ind+1) % 3

    rv = 1+np.random.uniform(-step_size, step_size)
    #print(rv)
    #rv = 1+0.5
    #transform = np.eye(3)
    new_cell[rnd_v1_ind,rnd_v1_ind] *= rv
    new_cell[rnd_v2_ind,rnd_v2_ind] *= (1/rv)
    
    delta_H = new_cell - cell
    print(delta_H)
    
    

    alk.alkane_change_box(int(ibox),delta_H, aspect_ratio_limit = 0.8)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    angle_limit_rad = angle_limit*np.pi/180
    
    if min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif alk.alkane_check_chain_overlap(int(ibox)):
        boltz = 0   
    else:
        boltz = 1
    
    return boltz, delta_H




        
def MC_run(ns_data, sweeps, move_ratio, ibox, volume_limit = sys.float_info.max, return_ase = False):    
    
        
    moves_accepted = np.zeros(6)
    moves_attempted = np.zeros(6)
    nbeads = alk.alkane_get_nbeads()
    nchains = alk.alkane_get_nchains()
    
    isweeps = 0
    pressure = 0
    move_prob = np.cumsum(move_ratio)/np.sum(move_ratio)
    
    if nbeads == 1:
        moves_per_sweep = nchains+7 
    elif nbeads <=3:
        moves_per_sweep = 2*nchains+7
    else:
        moves_per_sweep = (nbeads-1)*nchains+7 
    #moves_per_sweep selected such that every degree of freedom should 
    #be changed once on average when a sweep is performed.

        
    while isweeps < sweeps:
        imove=0
        while imove< moves_per_sweep:
#             clone_walker(ibox,nwalkers+2)#backup box
            ichain = np.random.randint(0, high=nchains) # picks a chain at random
            #should it be from 0 to nchains?
            current_chain = alk.alkane_get_chain(ichain+1, int(ibox))

            backup_chain = copy.deepcopy(current_chain)
            xi = np.random.random()
            if xi < move_prob[ivol]:
                # Attempt a volume move
                itype = ivol
                #clone_walker(ibox, volume_box) #backing up the volume
                boltz = alk.alkane_box_resize(pressure, int(ibox), 0)
                moves_attempted[itype] += 1
            elif xi < move_prob[itrans]:
                # Attempt a translation move
                itype = itrans
                boltz = alk.alkane_translate_chain(ichain+1, int(ibox))
                moves_attempted[itrans] += 1
            elif xi < move_prob[irot]:
                # Attempt a rotation move
                itype = irot
                boltz, quat = alk.alkane_rotate_chain(ichain+1, int(ibox), 0)
                moves_attempted[itype] += 1
            elif xi < move_prob[idih]:
                # Attempt a dihedral angle move
                itype = idih
                boltz, bead1, angle = alk.alkane_bond_rotate(ichain+1, int(ibox), 1)
                moves_attempted[itype] += 1
            elif xi < move_prob[ishear]:
                # Attempt a shear move
                itype = ishear
                boltz, delta_H = box_shear_step(ibox, ns_data)
                moves_attempted[itype] += 1
            else:
                # Attempt a stretch move
                itype = istr
                boltz, delta_H = box_shear_step(ibox, ns_data)
                moves_attempted[itype] += 1


            #Check which type of move and whether or not to accept
                    
            if (itype==ivol):
                new_volume = alk.box_compute_volume(int(ibox))
                if (np.random.random() < boltz) and (new_volume - volume_limit) < sys.float_info.epsilon:
                    moves_accepted[itype]+=1

                else:
                    #revert volume move
                    #clone_walker(volume_box, ibox)
                    dumboltz = alk.alkane_box_resize(pressure, int(ibox), 1)

            
            elif(itype == ishear or itype == istr):
                if boltz:
                    moves_accepted[itype]+=1

                
                else:
                    neg_delta_H = -1.0*delta_H
                    #print(neg_delta_H)
                    alk.alkane_change_box(int(ibox), neg_delta_H)

                

                
            else:
                if (np.random.random() < boltz):
                #accept the move
                    moves_accepted[itype]+=1


                else:
                    #reject the move
                    for ibead in range(nbeads):
                        current_chain[ibead] = backup_chain[ibead]

            imove += 1
        isweeps +=1
    moves_attempted = np.where(moves_attempted == 0, 1, moves_attempted)
    moves_acceptance_rate = moves_accepted/moves_attempted

    if return_ase:
        atoms =  mk_ase_config(ibox, nbeads, nchains, scaling=1)

        return atoms
    else:
        return alk.box_compute_volume(int(ibox)), moves_acceptance_rate

def clone_walker(ibox_source,ibox_clone):
    
    nbeads  = alk.alkane_get_nbeads()
    nchains = alk.alkane_get_nchains()
    



    cell = alk.box_get_cell(ibox_source)

    
    alk.box_set_cell(ibox_clone,cell)
    for ichain in range(1,nchains+1):
        original_chain = alk.alkane_get_chain(ichain,ibox_source)
        original_chain_copy = copy.deepcopy(original_chain)
        clone_chain = alk.alkane_get_chain(ichain,ibox_clone)
        for ibead in range(nbeads):
            clone_chain[ibead][:] = original_chain_copy[ibead][:]
    
    
    
    
def adjust_dv(ns_data,ibox,active_box, lower_bound,upper_bound, volume_limit, min_dv = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [1,0,0,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box, volume_limit)
    
    r = acceptance_rate[ivol]
    
    if r < lower_bound:
        alk.alkane_set_dv_max(max(alk.alkane_get_dv_max()/equil_factor, min_dv))
    elif r > upper_bound:
        alk.alkane_set_dv_max(min(alk.alkane_get_dv_max()*equil_factor, 10.0))
    return r



def adjust_dr(ns_data,ibox,active_box, lower_bound,upper_bound, min_dr = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,1,0,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[itrans]
    
    if r < lower_bound:
        alk.alkane_set_dr_max(max(alk.alkane_get_dr_max()/equil_factor,min_dr))
    elif r > upper_bound:
        alk.alkane_set_dr_max(min(alk.alkane_get_dr_max()*equil_factor,10.0))
    return r

def adjust_dt(ns_data,ibox,active_box, lower_bound,upper_bound, min_dt = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,1,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[irot]
    
    if r < lower_bound:
        alk.alkane_set_dt_max(max(alk.alkane_get_dt_max()/equil_factor, min_dt))
    elif r > upper_bound:
        alk.alkane_set_dt_max(alk.alkane_get_dt_max()*equil_factor)
    return r

def adjust_dh(ns_data,ibox,active_box, lower_bound,upper_bound, min_dh = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,1,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[idih]
    
    if r < lower_bound:
        alk.alkane_set_dh_max(max(alk.alkane_get_dh_max()/equil_factor,min_dh))
    elif r > upper_bound:
        alk.alkane_set_dh_max(alk.alkane_get_dh_max()*equil_factor)
    return r

def adjust_dshear(ns_data, ibox,active_box, lower_bound,upper_bound, min_dshear = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,1,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data, sweeps, move_ratio, active_box)
    
    r = acceptance_rate[ishear]
    
    
    if r < lower_bound:
        ns_data.shear_step_max = max((ns_data.shear_step_max*(1/equil_factor), min_dshear))
    elif r > upper_bound:
        ns_data.shear_step_max = ns_data.shear_step_max*equil_factor
    return r

def adjust_dstretch(ns_data, ibox,active_box, lower_bound,upper_bound, min_dstretch = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,0,1]
    clone_walker(ibox,active_box)
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[istr]
    
    
    if r < lower_bound:
        ns_data.stretch_step_max = max((ns_data.stretch_step_max*(1.0/equil_factor),min_dstretch))
    elif r > upper_bound:
        ns_data.stretch_step_max = ns_data.stretch_step_max*equil_factor
    return r

def perturb_initial_configs(ns_data, move_ratio, walk_length = 20):
    
    """ Runs a number of Monte Carlo steps on every simulation box, using the move_ratio assigned to it,
    Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume"""

    nwalkers = ns_data.parameters.nwalkers

    volumes = {}
    start_volumes = []
    for ibox in range(1,nwalkers+1):
        volumes[ibox], rate = MC_run(ns_data, walk_length, move_ratio, ibox)


    #overlap check
    overlap_check = np.zeros(nwalkers)
    for ibox in range(1,nwalkers+1):
        overlap_check[alk.alkane_check_chain_overlap(int(ibox))]
        

    return volumes


def celltoxmolstring(atoms):
    
    """Outputs the cell of an atoms object as a string in a format which can be used to write .xmol files
    
    Arguments:
        atoms: Atoms object for which a cell needs to be returned
        
    Returns:
        cellstring: String containing the three vectors which compose the cell of the atoms object"""
    
    cell = atoms.cell
    
    if cell.size == 3:
        cell *= np.eye(3)
    
    cellstring = np.array2string(atoms.cell.flatten(),
                                  max_line_width=100,
                                  formatter = {'float_kind':lambda x: "%.6f" % x})[1:-1]
    
    return cellstring

def write_configs_to_hdf(ns_data, current_iter, filename = None):
    """Write the coordinates of all hs_alkane boxes to a file. Default behavior is to overwrite previous data, 

    Arguments:
        filename (default = "Will use the filename provided by ns_data"): File to write atoms objects to.

    """
    if filename is None:
        filename = ns_data.restart_filename



    nwalkers = ns_data.parameters.nwalkers
    nbeads = ns_data.parameters.nbeads
    nchains = ns_data.parameters.nchains

    f = h5py.File(filename, "w")

    f.attrs.create("nbeads", nbeads)
    f.attrs.create("nchains", nchains)
    f.attrs.create("nwalkers", nwalkers)
    f.attrs.create("prev_iters", current_iter)
    f.attrs.create("sweeps", ns_data.parameters.walklength)
    

    for iwalker in range(1,nwalkers+1):
        groupname = f"walker_{iwalker:04d}"
        tempgrp = f.create_group(groupname)
        coords = tempgrp.create_dataset("coordinates",(nbeads*nchains,3),dtype="float64")
        unitcell = tempgrp.create_dataset("unitcell",(3,3),dtype="float64")

        unitcell[:] = alk.box_get_cell(iwalker)
        for ichain in range(nchains):
            chain = alk.alkane_get_chain(ichain+1,iwalker)
            coords[ichain*nbeads:ichain*nbeads+(nbeads), :] = chain

    f.close()

def adjust_mc_steps(ns_data, clone, active_box, volume_limit):

        # nbeads = ns_data.parameters.nbeads
        low_acc_rate = ns_data.low_acc_rate
        high_acc_rate = ns_data.high_acc_rate

        rates=[0,0,0,0,0,0]

        rates[0] = adjust_dv(ns_data,clone,active_box,low_acc_rate,high_acc_rate, volume_limit)
        rates[1] = adjust_dr(ns_data,clone,active_box,low_acc_rate,high_acc_rate)
        if ns_data.parameters.nbeads >= 2:
            rates[2] = adjust_dt(ns_data,clone,active_box,low_acc_rate,high_acc_rate)
        if ns_data.parameters.nbeads >= 4:
            rates[3] = adjust_dh(ns_data,clone,active_box,low_acc_rate,high_acc_rate)
        rates[4] = adjust_dshear(ns_data,clone,active_box,low_acc_rate,high_acc_rate)
        rates[5] = adjust_dstretch(ns_data,clone,active_box,low_acc_rate,high_acc_rate)


        return rates
    
def perform_ns_run(ns_data, iterations, prev_iters = 0, move_ratio = None, verbose = False):

    processes = ns_data.parameters.processes

    energies_file = open(ns_data.energies_filename, "a+")


    if verbose:
        print("Starting Nested Sampling run")
        print("iter  volume")
            
            
    active_box = ns_data.active_box
    # pool = mp.Pool(processes = processes)

    for i in range(prev_iters, prev_iters+iterations):


        index_max = ns_data.max_vol_index() #selecting the walker with highest volume
        volume_limit = ns_data.volumes[index_max] #setting volume limit


        energies_file.write(f"{i} {ns_data.volumes[index_max]} {ns_data.volumes[index_max]} \n") #writing this volume to output



        clone = np.random.randint(1,ns_data.parameters.nwalkers+1) #selecting which walker to clone

        clone_walker(clone, index_max)


        if i%ns_data.vis_interval == 0:
            ns_data.write_to_extxyz()
        
        if i%ns_data.mc_adjust_interval == 0:
            adjust_mc_steps(ns_data, clone, active_box, volume_limit)

        #Parallelism
        active_walkers = np.random.choice(np.setdiff1d(np.arange(1,ns_data.parameters.nwalkers+1),index_max,True), processes-1, False)
        active_walkers = np.append(clone, active_walkers)
        multiwalk = partial(MC_run, ns_data, ns_data.parameters.walklength, move_ratio, volume_limit = volume_limit, return_ase=True)



        # walk_output = map(MC_run, repeat(ns_data), repeat(ns_data.sweeps_per_walk),
        #                              repeat(move_ratio),active_walkers, repeat(volume_limit))

        pool = mp.Pool(processes)


        walk_output = pool.map(multiwalk,active_walkers)

        pool.close()
        pool.join()


        new_configs = [j for j in walk_output]


        # ns_data.volumes[active_walkers] = new_volumes
        # for j,walker in enumerate(active_walkers) :
        #     ns_data.volumes[walker] = new_volumes[j]

        for j, walker in enumerate(active_walkers):
            import_ase_to_ibox(new_configs[j], walker, ns_data)
            new_volume = alk.box_compute_volume(int(walker))
            ns_data.volumes[walker] = new_volume



        
        


        # new_config = MC_run(ns_data,ns_data.walklength, move_ratio, index_max,volume_limit, return_ase=True)

        # import_ase_to_ibox(new_config, index_max, ns_data)

        # new_volume=alk.box_compute_volume(index_max)

        
        # clone_walker(clone,active_box) #copies the ibox from first argument onto the second one.
        
        
        
        

        ns_data.volumes[index_max] = new_volume #updating the volume of the walker
        # clone_walker(active_box, index_max)
        if i%ns_data.print_interval == 0:
            #print(rates)
            print(f"{i:0>4}  {volume_limit}")


        if i%ns_data.restart_interval == 0:
            write_configs_to_hdf(ns_data,ns_data.restart_filename)
            if ns_data.time_remaining() < 1200:

                print("Out of allocated time, writing to file and exiting")

                energies_file.close()
                # pool.close()
                # pool.join()
                sys.exit()
    energies_file.close()
    # pool.close()
    # pool.join()

def import_ase_to_ibox(atoms, ibox, ns_data):
    """Inputs an ASE atoms object into a simulation cell."""


    cell_vectors = atoms.cell

    if cell_vectors.size == 3:
        cell_vectors *= np.eye(3)
    alk.box_set_cell(int(ibox),cell_vectors)

    positions = atoms.get_positions()

    for ichain in range(1,ns_data.parameters.nchains+1):
        chain = alk.alkane_get_chain(ichain,int(ibox))
        for ibead in range(ns_data.parameters.nbeads):
            chain[ibead][:] = positions[(ichain-1)*ns_data.parameters.nbeads+ibead][:]


    return


def initialise_sim_cells(ns_data):

    """Initialise hs_alkane cells"""
        
    alk.box_set_num_boxes(ns_data.parameters.nwalkers+1) #nwalkers+2 if debugging
    alk.box_initialise()
    alk.box_set_pbc(1)
    alk.alkane_set_nchains(ns_data.parameters.nchains) 
    alk.alkane_set_nbeads(ns_data.parameters.nbeads)    
    alk.alkane_initialise()           
    alk.box_set_isotropic(1)
    alk.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
    alk.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate thes





