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
from ase import io


#Setting constants 

#Take this out at some point, its bad
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']

ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5


class mc_step_sizes:

    def __init__(self):
        self.shear_step_max = float(0.5)
        self.stretch_step_max = float(0.5)


        alk.alkane_set_dr_max(0.65)
        alk.alkane_set_dt_max(0.43)
        alk.alkane_set_dh_max(0.4)
        alk.alkane_set_dv_max(0.5)

        self.vol_step_max = alk.alkane_get_dv_max()
        self.trans_step_max = alk.alkane_get_dr_max()
        self.rot_step_max = alk.alkane_get_dt_max()
        self.dih_step_max = alk.alkane_get_dh_max()
        

    def set_dshear(self,dx):
        self.shear_step_max = dx
        return
    
    def set_dstretch(self,dy):
        self.stretch_step_max = dy
        return

    def set_dv(self,dv):
        self.vol_step_max = dv
        alk.alkane_set_dv_max(dv)
        return

    def set_dr(self,dr):
        self.trans_step_max = dr
        alk.alkane_set_dr_max(dr)
        return

    def set_dt(self,dt):
        self.rot_step_max = dt
        alk.alkane_set_dt_max(dt)
        return

    def set_dh(self,dh):
        self.dih_step_max = dh
        alk.alkane_set_dh_max(dh)
        return
    
    def update_dv(self):
        alk.alkane_set_dv_max(self.vol_step_max)
        return

    def update_dr(self):
        alk.alkane_set_dr_max(self.trans_step_max)
        return
    
    def update_dt(self):
        alk.alkane_set_dt_max(self.rot_step_max)
        return

    def update_dh(self):
        alk.alkane_set_dh_max(self.dih_step_max)
        return

    def update_steps(self):
        self.update_dv()
        self.update_dt()
        self.update_dr()
        self.update_dh()

    
        





class ns_info:

    """Object containing most of the parameters required to perform the nested sampling simulation."""

    def __init__(self, parameters, from_restart = False):

        self.step_sizes=mc_step_sizes()
        self.parameters = parameters


        self.low_acc_rate = float(0.2)
        self.high_acc_rate = float(0.5)
        
        self.iter_ = int(0)

        self.from_restart = from_restart

        self.start_time = timer()

        self.active_box = parameters.nwalkers+1

    # # def __getattr__(self, item):
    #     tr
    # #     return getattr(self.parameters, item)
    
    
    def time_elapsed(self):
        """Returns the amount of time elapsed since this object was created:
        """
        current_time = timer()
        self.time_elapsed_ = current_time - self.start_time
        return self.time_elapsed_

    def time_remaining(self):
        """Returns how much time is remaining out of the alloted time for the simulation, assuming this object was created at the beginning of the program"""
        time_remaining_ = self.parameters.allotted_time - self.time_elapsed()

        return time_remaining_

    




    def set_dshear_max(self,dshear_max):
        """Set maximum size for shear moves:
            Arguments:
                dshear_max: The maximum size for shear moves to be set."""
        self.step_sizes.shear_step_max = dshear_max
        return

    def set_dstretch_max(self,dstretch_max):
        """Set maximum size for stretch moves:
            Arguments:
                dstretch_max: The maximum size for stretch moves to be set."""
        self.step_sizes.stretch_step_max = dstretch_max
        return

    def set_acc_rate_range(self,acc_range):
        """Set the lower and upper limit for the acceptance rates when adjusting the size of the MC moves.
            Arguments:
                acc_range: Array_like containing two floats, the first of which will be the lower limit for the move acceptance rate, and 
                           the second of which will be the upper limit for the acceptance rate."""
        if len(acc_range) != 2:
            raise IndexError("'acc_range' should have two values in it")
        if acc_range[0]>acc_range[1]:
            raise FloatingPointError("'acc_range[0]' should be less than 'acc_range[1]'")
        self.low_acc_rate = acc_range[0]
        self.high_acc_rate = acc_range[1]
        return

    def set_intervals(self,mc_adjust_interval = None,vis_interval = None, 
                        restart_interval = int(5e3), print_interval = 100):
        """Set the intervals at which various operations should be performed during the nested sampling run.
        Arguments:
            mc_adjust_interval: Set the interval at which the size of the Monte Carlo moves is adjusted. Defaults to nwalkers/2 rounded down.
            vis_interval: Set the interval at which a configuration is written to an extxyz file. Defaults to """
        if mc_adjust_interval is None:
            mc_adjust_interval = max(1,self.parameters.nwalkers//2)

        self.mc_adjust_interval = mc_adjust_interval
        
        if vis_interval is None:
            vis_interval = mc_adjust_interval

        self.vis_interval = vis_interval

        self.restart_interval = restart_interval

        self.print_interval = print_interval
        



    def check_overlaps(self):
        """Checks if overlaps are present in any of the simulation boxes.
        
            Returns:
                Overlaps: If 1, overlaps are present in one of the simulation boxes. If 0, no overlaps are present."""
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
        Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume.
        
            Arguments:
                move_ratio: Which move ratio to use during the Monte Carlo walk.
                walk_length: How many 'sweeps' to use in the Monte Carlo walk.
                
            Returns:
                self.volumes: A dictionary containing the volumes of all simulation boxes, with the keys corresponding the simulation box number, ibox."""

        nwalkers = self.parameters.nwalkers

        self.volumes = {}
        #start_volumes = []
        for ibox in range(1,nwalkers+1):
            self.volumes[ibox], _ = MC_run(self, walk_length, move_ratio, ibox)

        #overlap check
        self.check_overlaps()

        return self.volumes


    def set_directory(self,path="./"):
        """Set the path to working directory.
            Arguments:
                path: Directory where files will be written to."""

        self.pwd = path

        self.restart_filename = f"{path}restart.hdf5"
        self.traj_filename = f"{path}traj.extxyz"
        self.energies_filename = f"{path}energies.txt"

        return

        # self.energies_file = open(self.energies_filename, "a+")

        


            


    def load_volumes(self):
        """Construct a dictionary containing the volumes of all the simulation boxes in hs_alkane.
            Returns:
                self.volumes: A dictionary containing the volumes of all the simulation boxes."""
        self.volumes = {}
        for ibox in range(1,self.parameters.nwalkers+1):
            self.volumes[ibox] = alk.box_compute_volume(int(ibox))

        #overlap check
        self.check_overlaps()

        return self.volumes

    def write_to_extxyz(self, ibox=None):
        """Writes a single simulation box to file.
            Arguments:
                ibox: Simulation box to write. If none, the largest simulation box is used."""

        if ibox is None:
            ibox = self.max_vol_index()

        
        max_vol_config = mk_ase_config(ibox,self.parameters.nbeads,self.parameters.nchains)
        max_vol_config.wrap()

        io.write(self.traj_filename, max_vol_config, append = True)
        return

    def max_vol_index(self):
        """ Returns index of largest simulation box.
            Returns: Key corresponding to largest simulation box."""
        self.max_vol_index_ = max(self.volumes, key=self.volumes.get)
        return self.max_vol_index_

    def write_all_to_extxyz(self, filename = "dump.extxyz"):
        """ Writes all simulation boxes to an extxyz file:
        Arguments:
            filename: Name of file to output configurations to."""
        
        for i in range(1,self.parameters.nwalkers+1):
            atoms = mk_ase_config(i,self.parameters.nbeads,self.parameters.nchains)
            io.write(filename, atoms,format="extxyz", append = True)

        return


def mk_ase_config(ibox, Nbeads, Nchains, scaling = 3.75):
    """Uses the current state of the alkane model to construct an ASE atoms object.
        Arguments:
            ibox: Simulation box to convert to ASE object.
            Nbeads: Number of beads per chain in the configuration.
            Nchains: Number of chains in the configuration.
            scaling: A scaling factor to increase or decrease the size of the system by a constant. 3.75 is used as a. 
                    default as this produces configurations which approximate carbon as they have similar C-C distances.
        Returns:
            box_config: ASE atoms object"""
    
    # Create and populate ASE object
    model_positions = np.empty([Nchains*Nbeads, 3])
    cell_vectors = alk.box_get_cell(int(ibox))
   
    for ichain in range(0, Nchains):
        model_positions[Nbeads*ichain:Nbeads*ichain+Nbeads] = alk.alkane_get_chain(ichain+1, int(ibox))
    
    confstring = "C"+str(Nbeads*Nchains)
    
    box_config = Atoms(confstring, positions=model_positions*scaling, pbc=True, cell=cell_vectors*scaling)


    return box_config  # Returns ASE atom object

def vis_chains(vis_config, nbeads, nchains):
    """Takes an ASE atoms object or list thereof and creates a customised ngl viewer 
       with appropriate settings for our bulk alkane chain models.
       Arguments:
           vis_config: ASE object containing the configuration to be visualised.
           nbeads: Number of beads per chain in the configuration.
           nchains: Number of chains in the configuration.
        Returns:
            v: An ngl view instance, viewable in jupyter notebooks."""

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
    """Returns the shortest distance between two parallel faces, scaled such that the cell has a volume of 1.
    Arguments:
        ibox: Simulation box for which the min_aspect_ratio should be calculated for.
    Returns:
        min_aspect_ratio/np.cbrt(vol), A float representing the scaled shortest distance.
        
    """
    vol = alk.box_compute_volume(int(ibox))
    cell = alk.box_get_cell(int(ibox)).copy()
    #test_dist = []

    min_aspect_ratio = sys.float_info.max
    
    for i in range(3):
        vi = cell[i,:]
        vnorm_hat = np.cross(cell[(i+1)%3,:],cell[(i+2)%3,:])
        vnorm_hat = vnorm_hat/(np.sqrt(np.dot(vnorm_hat,vnorm_hat)))
        min_aspect_ratio = min(min_aspect_ratio, abs(np.dot(vnorm_hat,vi)))
        #test_dist.append(abs(np.dot(vnorm_hat,vi)))
        
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


def box_shear_step(ibox, ns_data, aspect_ratio_limit = 0.8, angle_limit = 45):
    """Perform a box shear move on a simulation box.
    Arguments:
        ibox: Simulation box on which to perform the box shear move.
        ns_data: ns_info object containing simulation parameter information.
        aspect_ratio_limit: Smallest allowed distance between parallel faces once normalised to unit volume, with larger values being more cubelike.
        angle_limit: Smallest allowed angle in degrees between two adjacent faces, to prevent the possibly squashing the unit cell.
    Returns:
        boltz: 0 if the proposed step has been rejected for being invalid, 1 if it is accepted.
        delta_H: Change in the unit cell, used in case the change in the cell should be reverted."""


    step_size = ns_data.step_sizes.shear_step_max
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
    
    alk.alkane_change_box(int(ibox),delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    
    angle_limit_rad = angle_limit*np.pi/180
    

    #reject due to poor aspect ratio
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

def box_stretch_step(ibox,ns_data, aspect_ratio_limit = 0.8, angle_limit = 45):    
    """Perform a box stretch move on a simulation box.
    Arguments:
        ibox: Simulation box on which to perform the box shear move.
        ns_data: ns_info object containing simulation parameter information.
        aspect_ratio_limit: Smallest allowed distance between parallel faces once normalised to unit volume, with larger values being more cubelike.
        angle_limit: Smallest allowed angle in degrees between two adjacent faces, to prevent the possibly squashing the unit cell.
    Returns:
        boltz: 0 if the proposed step has been rejected for being invalid, 1 if it is accepted.
        delta_H: Change in the unit cell, used in case the change in the cell should be reverted."""

    step_size = ns_data.step_sizes.stretch_step_max
    cell = alk.box_get_cell(int(ibox))
    new_cell = cell.copy()
    rnd_v1_ind = np.random.randint(0, 3)
    rnd_v2_ind = np.random.randint(0, 3)
    if rnd_v1_ind == rnd_v2_ind:
        rnd_v2_ind = (rnd_v2_ind+1) % 3

    rv = np.random.uniform(-step_size, step_size)
    #print(rv)
    #rv = 1+0.5
    #transform = np.eye(3)
    new_cell[rnd_v1_ind] *= np.exp(rv)
    new_cell[rnd_v2_ind] *= np.exp(-rv)
    
    delta_H = new_cell - cell
    
    

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

    ns_data.step_sizes.update_steps()

    # print(alk.alkane_get_dv_max(), alk.alkane_get_dr_max())
    
        
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
                boltz, delta_H = box_stretch_step(ibox, ns_data)
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
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box, volume_limit)
    
    r = acceptance_rate[ivol]
    
    if r < lower_bound:
        ns_data.step_sizes.set_dv(max(alk.alkane_get_dv_max()/equil_factor, min_dv))
    elif r > upper_bound:
        ns_data.step_sizes.set_dv(min(alk.alkane_get_dv_max()*equil_factor, 10.0))

    
    return r



def adjust_dr(ns_data,ibox,active_box, lower_bound,upper_bound, min_dr = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,1,0,0,0,0]
    clone_walker(ibox,active_box)
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[itrans]
    
    if r < lower_bound:
        ns_data.step_sizes.set_dr(max(alk.alkane_get_dr_max()/equil_factor,min_dr))
    elif r > upper_bound:
        ns_data.step_sizes.set_dr(min(alk.alkane_get_dr_max()*equil_factor,10.0))
    return r

def adjust_dt(ns_data,ibox,active_box, lower_bound,upper_bound, min_dt = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,1,0,0,0]
    clone_walker(ibox,active_box)
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[irot]
    
    if r < lower_bound:
        ns_data.step_sizes.set_dt(max(alk.alkane_get_dt_max()/equil_factor, min_dt))
    elif r > upper_bound:
        ns_data.step_sizes.set_dt(alk.alkane_get_dt_max()*equil_factor)
    return r

def adjust_dh(ns_data,ibox,active_box, lower_bound,upper_bound, min_dh = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,1,0,0]
    clone_walker(ibox,active_box)
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[idih]
    
    if r < lower_bound:
        ns_data.step_sizes.set_dh(max(alk.alkane_get_dh_max()/equil_factor,min_dh))
    elif r > upper_bound:
        ns_data.step_sizes.set_dh(alk.alkane_get_dh_max()*equil_factor)
    return r

def adjust_dshear(ns_data, ibox,active_box, lower_bound,upper_bound, min_dshear = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,1,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(ns_data, sweeps, move_ratio, active_box)
    
    r = acceptance_rate[ishear]
    
    
    if r < lower_bound:
        ns_data.step_sizes.set_dshear(max((ns_data.step_sizes.shear_step_max*(1/equil_factor), min_dshear)))
    elif r > upper_bound:
        ns_data.step_sizes.set_dshear(ns_data.step_sizes.shear_step_max*equil_factor)
    return r

def adjust_dstretch(ns_data, ibox,active_box, lower_bound,upper_bound, min_dstretch = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,0,1]
    clone_walker(ibox,active_box)
    _,acceptance_rate = MC_run(ns_data,sweeps, move_ratio, active_box)
    
    r = acceptance_rate[istr]
    
    
    if r < lower_bound:
        ns_data.step_sizes.set_dstretch(max((ns_data.step_sizes.stretch_step_max*(1.0/equil_factor),min_dstretch)))
    elif r > upper_bound:
        ns_data.step_sizes.set_dstretch(ns_data.step_sizes.stretch_step_max*equil_factor)
    return r

def perturb_initial_configs(ns_data, move_ratio, walk_length = 20):
    
    """ Runs a number of Monte Carlo steps on every simulation box, using the move_ratio assigned to it,
    Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume."""

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
    
    """Outputs the cell of an atoms object as a string in a format which can be used to write .xmol files.
    
    Arguments:
        atoms: Atoms object for which a cell needs to be returned.
        
    Returns:
        cellstring: String containing the three vectors which compose the cell of the atoms object."""
    
    cell = atoms.cell
    
    if cell.size == 3:
        cell *= np.eye(3)
    
    cellstring = np.array2string(atoms.cell.flatten(),
                                  max_line_width=100,
                                  formatter = {'float_kind':lambda x: "%.6f" % x})[1:-1]
    
    return cellstring

def write_configs_to_hdf(ns_data, current_iter, filename = None):
    """Write the coordinates of all hs_alkane boxes to a file. Default behavior is to overwrite previous data 

    Arguments:
        ns_data: Object containing simulation parameters
        current_iter: Which iteration of the nested sampling run is being written.
        filename: File to write atoms objects to. If none, use the filename provided by ns_data

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
    f.attrs.create("bondlength", ns_data.parameters.bondlength)
    f.attrs.create("bondangle", ns_data.parameters.bondangle)
    

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
        """Adjusts the size of the MC steps being performed on a box in order to correspond with a set acceptance rate, by performing MC runs on the boxes
        using only one move type.
        Arguments:
            ns_data: ns_info object containing simulation parameters.
            clone: Which simulation box to use to initialise the system on which the MC run is performed.
            active_box: Which simulation box to use to perform the runs on which stats are collected for adjusting the rate.
            volume_limit: The volume limit to be used when determining the acceptance rate of volume moves.
        Returns:
            rates: An array containing the acceptance rate for each type of MC move.

        """

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
    """Performs a nested sampling run using Monte Carlo walks, producing an energies.txt file, a restart file periodically and a trajectory file.
    Arguments:
        ns_data: ns_data object containing the simulation parameters
        iterations: Number of iterations to perform for the nested sampling run
        prev_iters: Number of previous iterations carried out during the nested sampling run. Use when working with a restart file.
        move_ratio: The ratio of moves to be used when performing Monte Carlo walks. Should be an array-like containing six values for each type of move,
                    corresponding with the moves: [volumes,'translation', rotation, dihedral, shear, stretch]
        verbose:    Whether or not to print data with runs (WIP)
                    """
        

    processes = ns_data.parameters.processes

    energies_file = open(ns_data.energies_filename, "a+")


    if verbose:
        print("Starting Nested Sampling run")
        print("iter  volume")
            
            
    active_box = ns_data.active_box
    pool = mp.Pool(processes = processes)

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

        active_walkers = np.random.choice(np.setdiff1d(np.arange(1,ns_data.parameters.nwalkers+1),index_max,True), processes-1, False)
        active_walkers = np.append(index_max, active_walkers)
        multiwalk = partial(ase_MC_run, ns_data=ns_data, sweeps=ns_data.parameters.walklength, move_ratio=move_ratio, volume_limit = volume_limit, return_ase=True,
                        ibox=active_box)


        walk_input = [mk_ase_config(j, ns_data.parameters.nbeads, ns_data.parameters.nchains, scaling=1) for j in active_walkers]
        #list of atomic configurations to walk


        walk_output = pool.map(multiwalk,walk_input)


        new_configs = [j for j in walk_output]

        for j, walker in enumerate(active_walkers):
            import_ase_to_ibox(new_configs[j], walker, ns_data)
            new_volume = alk.box_compute_volume(int(walker))
            ns_data.volumes[walker] = new_volume

        # clone_walker(active_box, index_max)
        if i%ns_data.print_interval == 0:
        # if i %10 ==0:
            #print(rates)
            print(f"{i:0>4}  {volume_limit}")

            
        if ns_data.time_remaining() < 0:

            print("Out of allocated time, writing to file and exiting")
            write_configs_to_hdf(ns_data,i,filename=ns_data.restart_filename)

            # energies_file.close()
            # pool.close()
            # pool.join()
            break
    energies_file.close()
    pool.close()
    pool.join()

def import_ase_to_ibox(atoms, ibox, ns_data):
    """Inputs an ASE atoms object into a simulation cell.
    Arguments:
        atoms: ns_data object containing the parameters for the simulation.
        ibox: Which `hs_alkane` simulation box to import the ase data into
        ns_data: ns_data object containing the simulation parameters
        """


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

    """Initialise hs_alkane cells
    
    Arguments:
        ns_data: ns_data object containing the parameters for the simulation."""
        
    alk.box_set_num_boxes(ns_data.parameters.nwalkers+1) #nwalkers+2 if debugging
    alk.box_initialise()
    alk.box_set_pbc(1)
    alk.alkane_set_nchains(ns_data.parameters.nchains)
    alk.alkane_set_nbeads(ns_data.parameters.nbeads)
    alk.alkane_initialise()           
    alk.box_set_isotropic(1)
    alk.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
    alk.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate these
    alk.alkane_set_bondlength(ns_data.parameters.bondlength)
    alk.alkane_set_bondangle(ns_data.parameters.bondangle)

def ase_MC_run(atoms, **kwargs):

    import_ase_to_ibox(atoms, kwargs["ibox"], kwargs["ns_data"])

    new_atoms = MC_run(**kwargs)

    return new_atoms






