# initial params/setting up empty boxes
import hs_alkane.alkane as alk
import numpy as np
import sys
import copy
from ase import io
from timeit import default_timer as timer
#import cProfile
from mpi4py import MPI

#for reading/writing



#for converting hs_alkane boxes to ASE atoms objects
from ase import Atoms
from ase.visualize import view
from ase import io


#Setting constants 

#Take this out at some point, its bad
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']

ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5


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

def box_shear_step(ibox, step_size, aspect_ratio_limit = 0.8, angle_limit = 60):
    """Perform a box shear move on a simulation box.
    Arguments:
        ibox: Simulation box on which to perform the box shear move.
        ns_data: ns_info object containing simulation parameter information.
        aspect_ratio_limit: Smallest allowed distance between parallel faces once normalised to unit volume, with larger values being more cubelike.
        angle_limit: Smallest allowed angle in degrees between two adjacent faces, to prevent the possibly squashing the unit cell.
    Returns:
        boltz: 0 if the proposed step has been rejected for being invalid, 1 if it is accepted.
        delta_H: Change in the unit cell, used in case the change in the cell should be reverted."""

    # pick random vector for shear direction
    #np.random.seed(10)
    rnd_vec_ind = int(np.floor(alk.random_uniform_random()*3))
    # turn other two into orthonormal pair
    #should I have a pair or should I have a single vector.
    other_vec_ind = list(range(3))
    other_vec_ind.remove(rnd_vec_ind)
    # orig_cell_copy = copy.deepcopy(alk.box_get_cell(int(ibox)))
    orig_cell_copy = alk.box_get_cell(int(ibox)).copy()

    
    v1 = orig_cell_copy[other_vec_ind[0],:]
    v2 = orig_cell_copy[other_vec_ind[1],:]  

    
    v1 /= np.sqrt(np.dot(v1,v1))
    v2 -= v1*np.dot(v1,v2) 
    v2 /= np.sqrt(np.dot(v2,v2))

    if np.isnan(np.sum(v1)) or np.isnan(np.sum(v2)):
        print(v1,v2)
        print(orig_cell_copy)
        sys.exit()



    
    # pick random magnitudes
    rv1 = np.random.uniform(-step_size, step_size)
    rv2 = np.random.uniform(-step_size, step_size)
#     rv1 = 2
#     rv2 = 0

    # create new cell and transformation matrix (matrix is additive)
    new_cell = orig_cell_copy.copy()

    new_cell[rnd_vec_ind,:] += rv1*v1 + rv2*v2

    delta_H = new_cell - orig_cell_copy
    
    alk.alkane_change_box(int(ibox),delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    
    angle_limit_rad = angle_limit*np.pi/180
    

    #reject due to poor aspect ratio
    if alk.box_min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif alk.alkane_check_chain_overlap(int(ibox)):
        boltz = 0   
    else:
        boltz = 1
    
    
    #bake rejection  due to shape in here as it becomes easier to fit with the rest of the code
    
    
    return boltz, delta_H

def box_stretch_step(ibox,step_size, aspect_ratio_limit = 0.8, angle_limit = 60):    
    """Perform a box stretch move on a simulation box.
    Arguments:
        ibox: Simulation box on which to perform the box shear move.
        ns_data: ns_info object containing simulation parameter information.
        aspect_ratio_limit: Smallest allowed distance between parallel faces once normalised to unit volume, with larger values being more cubelike.
        angle_limit: Smallest allowed angle in degrees between two adjacent faces, to prevent the possibly squashing the unit cell.
    Returns:
        boltz: 0 if the proposed step has been rejected for being invalid, 1 if it is accepted.
        delta_H: Change in the unit cell, used in case the change in the cell should be reverted."""

    cell = alk.box_get_cell(int(ibox))
    new_cell = cell.copy()
    rnd_v1_ind = int(np.floor(alk.random_uniform_random()*3))
    rnd_v2_ind = int(np.floor(alk.random_uniform_random()*3))
    if rnd_v1_ind == rnd_v2_ind:
        rnd_v2_ind = (rnd_v2_ind+1) % 3

    rv = np.random.uniform(-step_size, step_size)
    #print(rv)
    #rv = 1+0.5
    #transform = np.eye(3)
    new_cell[rnd_v1_ind] *= np.exp(rv)
    new_cell[rnd_v2_ind] *= np.exp(-rv)
    
    delta_H = new_cell - cell
    
    

    alk.alkane_change_box(int(ibox),delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    angle_limit_rad = angle_limit*np.pi/180
    
    if alk.box_min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif alk.alkane_check_chain_overlap(int(ibox)):
        boltz = 0   
    else:
        boltz = 1
    
    return boltz, delta_H

def MC_run(ns_data, sweeps, move_ratio, ibox, volume_limit = sys.float_info.max, return_ase = False,
dshear = 1.0, dstretch = 1.0, min_ang = 60, min_ar = 0.8):

    #ns_data.step_sizes.update_steps()

    # print(alk.alkane_get_dv_max(), alk.alkane_get_dr_max())
    moves_accepted = np.zeros(6)
    moves_attempted = np.zeros(6)
    nbeads = ns_data["nbeads"]
    nchains = ns_data["nchains"]
    
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
            ichain = int(np.floor(alk.random_uniform_random()*nchains)) # picks a chain at random
            #should it be from 0 to nchains?
            current_chain = alk.alkane_get_chain(ichain+1, int(ibox))

            backup_chain = current_chain.copy()
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
                boltz, delta_H = box_shear_step(ibox, dshear, min_ar, min_ang)
                moves_attempted[itype] += 1
            else:
                # Attempt a stretch move
                itype = istr
                boltz, delta_H = box_stretch_step(ibox, dstretch, min_ar,min_ang)
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
        clone_chain = alk.alkane_get_chain(ichain,ibox_clone)
        for ibead in range(nbeads):
            clone_chain[ibead][:] = original_chain[ibead][:]
      

def perturb_initial_configs(ns_data, move_ratio, walk_length = 20):
    
    """ Runs a number of Monte Carlo steps on every simulation box, using the move_ratio assigned to it,
    Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume."""

    nwalkers = ns_data["nwalkers"]

    volumes = {}
    start_volumes = []
    for ibox in range(1,nwalkers+1):
        volumes[ibox], rate = MC_run(ns_data, walk_length, move_ratio, ibox)


    #overlap check
    overlap_check = np.zeros(nwalkers)
    for ibox in range(1,nwalkers+1):
        overlap_check[alk.alkane_check_chain_overlap(int(ibox))]
        

    return volumes



def import_ase_to_ibox(atoms, ibox, ns_data, scaling = 1.0):
    """Inputs an ASE atoms object into a simulation cell.
    Arguments:
        atoms: ns_data object containing the parameters for the simulation.
        ibox: Which `hs_alkane` simulation box to import the ase data into
        ns_data: ns_data object containing the simulation parameters
        """

    try:
        nbeads = ns_data.parameters.nbeads
        nchains = ns_data.parameters.nchains
    except:
        nbeads = ns_data["nbeads"]
        nchains = ns_data["nchains"]


    cell_vectors = atoms.cell

    if cell_vectors.size == 3:
        cell_vectors *= np.eye(3)
    alk.box_set_cell(int(ibox),cell_vectors*scaling)

    positions = atoms.get_positions()

    for ichain in range(1,nchains+1):
        chain = alk.alkane_get_chain(ichain,int(ibox))
        for ibead in range(nbeads):
            chain[ibead][:] = positions[(ichain-1)*nbeads+ibead][:]*scaling


    return

def initialise_sim_cells(args, quiet):

    """Initialise hs_alkane cells
    
    Arguments:
        ns_data: ns_data object containing the parameters for the simulation."""


    # alk.random_set_random_seed(1)
    alk.box_set_quiet(quiet)
    alk.box_set_num_boxes(int(args["nwalkers"])) #nwalkers+2 if debugging
    alk.box_initialise()
    alk.box_set_pbc(1)
    alk.alkane_set_nchains(int(args["nchains"]))
    alk.alkane_set_nbeads(int(args["nbeads"]))
    alk.alkane_initialise()           
    alk.box_set_isotropic(1)
    alk.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
    alk.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate these
    alk.alkane_set_bondlength(float(args["bondlength"]))
    alk.alkane_set_bondangle(float(args["bondangle"]))



def default_move_ratio(ns_data):
    if "move_ratio" in ns_data:
        move_ratio = [float(i) for i in ns_data["move_ratio"].split(',')]
        assert len(move_ratio) == 6, "Move ratio array should be of length 6."
        pass
        return move_ratio
    move_ratio = np.zeros(6)
    move_ratio[ivol] = 1
    move_ratio[itrans] = 3.0*ns_data["nchains"]
    move_ratio[irot] = (2.0*ns_data["nchains"]) if ns_data["nbeads"] >= 2 else 0
    move_ratio[idih] = 1.0*max(((ns_data["nbeads"]-3.0)*(ns_data["nchains"]),0))
    move_ratio[ishear] = 3
    move_ratio[istr] = 3
    return move_ratio

def create_initial_configs(args, max_vol_per_atom = 15):
    cell_matrix = 0.999*np.eye(3)*np.cbrt(args["nbeads"]*args["nchains"]*max_vol_per_atom)#*np.random.uniform(0,1)
    for ibox in range(1,args["nwalkers"]+1):
        alk.box_set_cell(int(ibox),cell_matrix)
    populate_boxes(args)

def populate_boxes(args):
    ncopy = args["nchains"]
    for ibox in range(1,args["nwalkers"]+1):
        for ichain in range(1,ncopy+1):
            rb_factor = 0
            alk.alkane_set_nchains(ichain)
            overlap_flag = 1
            while rb_factor == 0:
                rb_factor, ifail = alk.alkane_grow_chain(ichain,int(ibox),1) 
                if ifail != 0:
                    rb_factor = 0

def write_all_to_extxyz(args,filename = "dump.extxyz"):
    """ Writes all simulation boxes to an extxyz file:
    Arguments:
        filename: Name of file to output configurations to."""
    
    for i in range(1,args["nwalkers"]+1):
        atoms = mk_ase_config(i,args["nbeads"],args["nchains"])
        io.write(filename, atoms,format="extxyz", append = True)

    return

def adjust_mc_steps(args,comm,move_ratio,vol_max,walklength = 10, lower_bound = 0.2, upper_bound=0.5, 
                      min_dstep=1e-5*np.ones(6), dv_max=10.0,dr_max=10.0, dshear = 1.0, dstretch=1.0): 

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
  

    size = comm.Get_size()
    rate = np.zeros(6)
    avg_rate = np.zeros_like(rate)
    mc_box = np.random.randint(args["nwalkers"])

    for i in range(6):
        move_ratio_matrix = np.eye(6)
        backup = mk_ase_config(mc_box+1,args["nbeads"],args["nchains"],scaling=1)
        if move_ratio[i] != 0:
            rate += MC_run(args,walklength, move_ratio_matrix[i],mc_box+1,vol_max,dshear=dshear, dstretch=dstretch)[1]
            import_ase_to_ibox(backup,mc_box+1,args)
    comm.Allreduce(rate,avg_rate,op=MPI.SUM)
    avg_rate = avg_rate/size
    if move_ratio[0] != 0:
        if avg_rate[0] < lower_bound:
            alk.alkane_set_dv_max(max(0.5*alk.alkane_get_dv_max(),min_dstep[0]))
        elif avg_rate[0] > upper_bound:
            alk.alkane_set_dv_max(min(2.0*alk.alkane_get_dv_max(),dv_max))
    if move_ratio[1] != 0:
        if avg_rate[1] < lower_bound:
            alk.alkane_set_dr_max(max(0.5*alk.alkane_get_dr_max(),min_dstep[1]))
        elif avg_rate[1] > upper_bound:
            alk.alkane_set_dr_max(min(2.0*alk.alkane_get_dr_max(),dr_max))
    if move_ratio[2] != 0:
        if avg_rate[2] < lower_bound:
            alk.alkane_set_dt_max(max(0.5*alk.alkane_get_dt_max(),min_dstep[2]))
        elif avg_rate[2] > upper_bound:
            alk.alkane_set_dt_max(2.0*alk.alkane_get_dt_max())
    if move_ratio[3] != 0:
        if avg_rate[3] < lower_bound:
            alk.alkane_set_dh_max(max(0.5*alk.alkane_get_dh_max(),min_dstep[3]))
        elif avg_rate[3] > upper_bound:
            alk.alkane_set_dh_max(2.0*alk.alkane_get_dh_max())
    if move_ratio[4] != 0:
        if avg_rate[4] < lower_bound:
            dshear = max(0.5*dshear,min_dstep[4])
        elif avg_rate[4] > upper_bound:
            dshear = 2.0*dshear
    if move_ratio[5] != 0:
        if avg_rate[5] < lower_bound:
            dstretch = max(0.5*dstretch,min_dstep[5])
        elif avg_rate[5] > upper_bound:
            dstretch  = 2.0*dstretch
    return avg_rate, dshear, dstretch

def signal_handler(signal, frame):
    global interrupted
    interrupted = True
