# initial params/setting up empty boxes
import hs_alkane.alkane as mdl
import numpy as np
import sys
import os
import copy

import argparse #parse arguments


#for converting hs_alkane boxes to ASE atoms objects
from ase import Atoms
#from ase.visualize import view
from ase.io import write as asewrite



#Setting constants
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']

ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5

shear_step_max = 0.5
stretch_step_max = 0.5



def mk_ase_config(ibox, Nbeads, Nchains, scaling = 3.75):
    'Uses the current state of the alkane model to construct an ASE atoms object'
    
    # Create and populate ASE object
    model_positions = np.empty([Nchains*Nbeads, 3])
    cell_vectors = mdl.box_get_cell(ibox)
   
    for ichain in range(0, Nchains):
        model_positions[Nbeads*ichain:Nbeads*ichain+Nbeads] = mdl.alkane_get_chain(ichain+1, ibox)
    
    confstring = "C"+str(Nbeads*Nchains)
    
    box_config = Atoms(confstring, positions=model_positions*scaling, pbc=True, cell=cell_vectors*scaling)


    return box_config  # Returns ASE atom object

def vis_chains(vis_config):
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
    vol = mdl.box_compute_volume(ibox)
    cell = mdl.box_get_cell(ibox).copy()
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
    cell = mdl.box_get_cell(ibox).copy()
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


def box_shear_step(ibox, step_size, shear_step_max = 0.5, aspect_ratio_limit = 0.8, angle_limit = 55):
    # pick random vector for shear direction
    #np.random.seed(10)
    rnd_vec_ind = np.random.randint(0, 3)
    # turn other two into orthonormal pair
    #should I have a pair or should I have a single vector.
    other_vec_ind = list(range(3))
    other_vec_ind.remove(rnd_vec_ind)
    orig_cell = copy.deepcopy(mdl.box_get_cell(ibox))
    orig_cell_copy = copy.deepcopy(orig_cell)

    
    v1 = orig_cell_copy[other_vec_ind[0],:]
    v2 = orig_cell_copy[other_vec_ind[1],:]
    
    v1 /= np.sqrt(np.dot(v1,v1))
    v2 -= v1*np.dot(v1,v2) 
    v2 /= np.sqrt(np.dot(v2,v2))

    
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

    

    
    mdl.alkane_change_box(ibox,delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    
    angle_limit_rad = angle_limit*np.pi/180
    
    if min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif mdl.alkane_check_chain_overlap(ibox):
        boltz = 0   
    else:
        boltz = 1
    
    
    #bake rejection  due to shape in here as it becomes easier to fit with the rest of the code
    
    
    return boltz, delta_H

def box_stretch_step(ibox, step_size, stretch_step_max = 0.5, aspect_ratio_limit = 0.8, angle_limit = 55):
    cell = mdl.box_get_cell(ibox)
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
    

    
    mdl.alkane_change_box(ibox,delta_H, aspect_ratio_limit = 0.8)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    
    angle_limit_rad = angle_limit*np.pi/180
    
    if min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif mdl.alkane_check_chain_overlap(ibox):
        boltz = 0   
    else:
        boltz = 1
    
    return boltz, delta_H




        
def MC_run(sweeps, move_ratio, ibox, volume_limit = sys.float_info.max):
    
    
        
    moves_accepted = np.zeros(6)
    moves_attempted = np.zeros(6)
    nbeads = mdl.alkane_get_nbeads()
    nchains = mdl.alkane_get_nchains()
    
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
            current_chain = mdl.alkane_get_chain(ichain+1, ibox)
            backup_chain = copy.deepcopy(current_chain)
            xi = np.random.random()
            if xi < move_prob[ivol]:
                # Attempt a volume move
                itype = ivol
                #clone_walker(ibox, volume_box) #backing up the volume
                boltz = mdl.alkane_box_resize(pressure, ibox, 0)
                moves_attempted[itype] += 1
            elif xi < move_prob[itrans]:
                # Attempt a translation move
                itype = itrans
                boltz = mdl.alkane_translate_chain(ichain+1, ibox)
                moves_attempted[itrans] += 1
            elif xi < move_prob[irot]:
                # Attempt a rotation move
                itype = irot
                boltz, quat = mdl.alkane_rotate_chain(ichain+1, ibox, 0)
                moves_attempted[itype] += 1
            elif xi < move_prob[idih]:
                # Attempt a dihedral angle move
                itype = idih
                boltz, bead1, angle = mdl.alkane_bond_rotate(ichain+1, ibox, 1)
                moves_attempted[itype] += 1
            elif xi < move_prob[ishear]:
                # Attempt a shear move
                itype = ishear
                boltz, delta_H = box_shear_step(ibox, shear_step_max)
                moves_attempted[itype] += 1
            else:
                # Attempt a stretch move
                itype = istr
                boltz, delta_H = box_shear_step(ibox, shear_step_max)
                moves_attempted[itype] += 1

#             if (np.random.random() < boltz):
#                 #accept the move
#                 moves_accepted[itype]+=1

#             else:
                
#                 # Reject move
#                 if (itype!=ivol):
#                     # Restore old chain if single chain move
#                     for ibead in range(nbeads):
#                         current_chain[ibead] = backup_chain[ibead]
#                 else:
#                     # Reset the box change - special fucntion for this.
#                     dumboltz = mdl.alkane_box_resize(pressure, ibox, 1)
                    
            if (itype==ivol):
                new_volume = mdl.box_compute_volume(ibox)
                if (np.random.random() < boltz) and (new_volume - volume_limit) < sys.float_info.epsilon:
                    moves_accepted[itype]+=1
#                     if nan_check(ibox):
#                         print(f"NAN accepted as {move_types[itype]} move")
#                         sys.exit()
                else:
                    #revert volume move
                    #clone_walker(volume_box, ibox)
                    dumboltz = mdl.alkane_box_resize(pressure, ibox, 1)
#                     if nan_check(ibox):
#                         print(f"NAN rejected trying to revert {move_types[itype]} move")
#                         sys.exit()
            
            elif(itype == ishear or itype == istr):
                if boltz:
                    moves_accepted[itype]+=1
#                     if nan_check(ibox):
#                         print(f"NAN accepted as {move_types[itype]} move")
#                         sys.exit()
                
                else:
                    neg_delta_H = -1.0*delta_H
                    #print(neg_delta_H)
                    mdl.alkane_change_box(ibox, neg_delta_H)
#                     if nan_check(ibox):
#                         print(f"NAN rejected trying to revert {move_types[itype]} move")
#                         sys.exit()
                

                
            else:
                if (np.random.random() < boltz):
                #accept the move
                    moves_accepted[itype]+=1
#                     if nan_check(ibox):
#                         print(f"NAN accepted as {move_types[itype]} move")
#                         print("current_chain", "backup chain")
#                         print(current_chain, backup_chain)
#                         sys.exit()

                else:
                    #reject the move
                    for ibead in range(nbeads):
                        current_chain[ibead] = backup_chain[ibead]
#                     if nan_check(ibox):
#                         print(f"NAN rejected trying to revert {move_types[itype]} move")
#                         print("current_chain", "backup chain")
#                         print(current_chain, backup_chain)
#                         sys.exit()
            imove += 1
        isweeps +=1
    moves_attempted = np.where(moves_attempted == 0, 1, moves_attempted)
    moves_acceptance_rate = moves_accepted/moves_attempted

    return mdl.box_compute_volume(ibox), moves_acceptance_rate

def clone_walker(ibox_original,ibox_clone):
    
    nbeads  = mdl.alkane_get_nbeads()
    nchains = mdl.alkane_get_nchains()
    
    cell = mdl.box_get_cell(ibox_original)
    
    mdl.box_set_cell(ibox_clone,cell)
    for ichain in range(1,nchains+1):
        original_chain = mdl.alkane_get_chain(ichain,ibox_original)
        original_chain_copy = copy.deepcopy(original_chain)
        clone_chain = mdl.alkane_get_chain(ichain,ibox_clone)
        for i in range(nbeads):
#             if np.isnan(original_chain_copy).any():
#                 print("original","clone")
#                 print(ibox_original, ibox_clone)
                
#                 print("clone_chain")
#                 print(clone_chain)
                
#                 print("original_chain_copy")
#                 print(original_chain_copy)
                
#                 print("original_chain")
#                 print(original_chain)
                
#                 sys.exit()
            clone_chain[i][:] = original_chain_copy[i][:]
    
    
    
    
def adjust_dv(ibox,active_box, lower_bound,upper_bound, volume_limit, min_dv = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [1,0,0,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box, volume_limit)
    
    r = acceptance_rate[ivol]
    
    if r < lower_bound:
        mdl.alkane_set_dv_max(max(mdl.alkane_get_dv_max()/equil_factor, min_dv))
    elif r > upper_bound:
        mdl.alkane_set_dv_max(min(mdl.alkane_get_dv_max()*equil_factor, 10.0))
    return r



def adjust_dr(ibox,active_box, lower_bound,upper_bound, min_dr = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,1,0,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[itrans]
    
    if r < lower_bound:
        mdl.alkane_set_dr_max(max(mdl.alkane_get_dr_max()/equil_factor,min_dr))
    elif r > upper_bound:
        mdl.alkane_set_dr_max(min(mdl.alkane_get_dr_max()*equil_factor,10.0))
    return r

def adjust_dt(ibox,active_box, lower_bound,upper_bound, min_dt = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,1,0,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[irot]
    
    if r < lower_bound:
        mdl.alkane_set_dt_max(max(mdl.alkane_get_dt_max()/equil_factor, min_dt))
    elif r > upper_bound:
        mdl.alkane_set_dt_max(mdl.alkane_get_dt_max()*equil_factor)
    return r

def adjust_dh(ibox,active_box, lower_bound,upper_bound, min_dh = 1e-10):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,1,0,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[idih]
    
    if r < lower_bound:
        mdl.alkane_set_dh_max(max(mdl.alkane_get_dh_max()/equil_factor,min_dh))
    elif r > upper_bound:
        mdl.alkane_set_dh_max(mdl.alkane_get_dh_max()*equil_factor)
    return r

def adjust_dshear(ibox,active_box, lower_bound,upper_bound, min_dshear = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,1,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[ishear]
    
    global shear_step_max
    
    if r < lower_bound:
        shear_step_max = max((shear_step_max*(1/equil_factor), min_dshear))
    elif r > upper_bound:
        shear_step_max = shear_step_max*equil_factor
    return r

def adjust_dstretch(ibox,active_box, lower_bound,upper_bound, min_dstretch = 1e-5):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,0,1]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[istr]
    
    global stretch_step_max
    
    if r < lower_bound:
        stretch_step_max = max((stretch_step_max*(1.0/equil_factor),min_dstretch))
    elif r > upper_bound:
        stretch_step_max = stretch_step_max*equil_factor
    return r

def populate_boxes(nwalkers,nchains):
    ncopy = nchains
    for ibox in range(1,nwalkers+1):
        for ichain in range(1,ncopy+1):
            rb_factor = 0
            mdl.alkane_set_nchains(ichain)
            overlap_flag = 1
            while rb_factor == 0:
                rb_factor, ifail = mdl.alkane_grow_chain(ichain,ibox,1)         
                if ifail != 0:
                    rb_factor = 0

def shake_initial_configs(nwalkers, move_ratio, walk_length = 20):
    
    """ Runs a number of Monte Carlo steps on every simulation box, using the move_ratio assigned to it,
    Checks for overlaps, and returns a dictionary which uses the number for each simulation box as the key for its volume"""


    volumes = {}
    start_volumes = []
    for ibox in range(1,nwalkers+1):
        volumes[ibox], rate = MC_run(walk_length, move_ratio, ibox)


    #overlap check
    overlap_check = np.zeros(nwalkers)
    for ibox in range(1,nwalkers+1):
        overlap_check[mdl.alkane_check_chain_overlap(ibox)]
        

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
