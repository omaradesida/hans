# initial params/setting up empty boxes
import hs_alkane.alkane as mdl
import numpy as np
import sys
import os
import copy


#visualise data
import matplotlib.pyplot as plt

#for converting hs_alkane boxes to ASE atoms objects
from ase import Atoms
#from ase.visualize import view
import ase

def mk_ase_config(ibox, Nbeads, Nchains):
    'Uses the current state of the alkane model to construct an ASE atoms object'
    
    # Create and populate ASE object
    model_positions = np.empty([Nchains*Nbeads, 3])
    cell_vectors = mdl.box_get_cell(ibox)
   
    for ichain in range(0, Nchains):
        model_positions[Nbeads*ichain:Nbeads*ichain+Nbeads] = mdl.alkane_get_chain(ichain+1, ibox)
    
    confstring = "C"+str(Nbeads*Nchains)
    
    box_config = Atoms(confstring, positions=model_positions*(1.5/0.4), pbc=True, cell=cell_vectors*(1.5/0.4))


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


def box_shear_step(ibox, step_size):
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

def box_stretch_step(ibox, step_size):
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
    

    
    mdl.alkane_change_box(ibox,delta_H)
    #transform = np.dot(np.linalg.inv(orig_cell), new_cell)
    if min_aspect_ratio(ibox) < aspect_ratio_limit:
        boltz = 0
    elif min_angle(ibox) < angle_limit_rad:
        boltz = 0
    elif mdl.alkane_check_chain_overlap(ibox):
        boltz = 0   
    else:
        boltz = 1
    
    return boltz, delta_H


#populating boxes
ncopy = nchains
for ibox in range(1,nwalkers+1):
    for ichain in range(1,ncopy+1):
        rb_factor = 0
        mdl.alkane_set_nchains(ichain)
        overlap_flag = 1
        while overlap_flag != 0:
            rb_factor, ifail = mdl.alkane_grow_chain(ichain,ibox,1)         
#             if ifail != 0:
#                 rb_factor = 0
            overlap_flag = mdl.alkane_check_chain_overlap(ibox)

        
def MC_run(sweeps, move_ratio, ibox, volume_limit = sys.float_info.max):
        
    moves_accepted = np.zeros(6)
    moves_attempted = np.zeros(6)
    nbeads = mdl.alkane_get_nbeads()
    nchains = mdl.alkane_get_nchains()
    isweeps = 0
    pressure = 0
    move_prob = np.cumsum(move_ratio)/np.sum(move_ratio)
    while isweeps < sweeps:
        imove=0
        while imove< nchains+7:
#             clone_walker(ibox,nwalkers+2)#backup box
            ichain = np.random.randint(0, high=nchains) # picks a chain at random
            #should it be from 0 to nchains?
            current_chain = mdl.alkane_get_chain(ichain+1, ibox)
            backup_chain = copy.deepcopy(current_chain)
            xi = np.random.random()
            if xi < move_prob[ivol]:
                # Attempt a volume move
                itype = ivol
                clone_walker(ibox, volume_box) #backing up the volume
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
                    clone_walker(volume_box, ibox)
                    #dumboltz = mdl.alkane_box_resize(pressure, ibox, 1)
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
    
    
    
    
def adjust_dv(ibox,active_box, lower_bound,upper_bound, volume_limit):
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



def adjust_dr(ibox,active_box, lower_bound,upper_bound):
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

def adjust_dt(ibox,active_box, lower_bound,upper_bound):
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

def adjust_dh(ibox,active_box, lower_bound,upper_bound):
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

def adjust_dshear(ibox,active_box, lower_bound,upper_bound):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,1,0]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[ishear]
    
    global shear_step_max
    
    if r < lower_bound:
        shear_step_max = max(shear_step_max*(1/equil_factor), min_dshear)
    elif r > upper_bound:
        shear_step_max = shear_step_max*equil_factor
    return r

def adjust_dstretch(ibox,active_box, lower_bound,upper_bound):
    equil_factor = 2
    sweeps = 20
    move_ratio = [0,0,0,0,0,1]
    clone_walker(ibox,active_box)
    vol,acceptance_rate = MC_run(sweeps, move_ratio, active_box)
    
    r = acceptance_rate[istr]
    
    global stretch_step_max
    
    if r < lower_bound:
        stretch_step_max = max(stretch_step_max*(1.0/equil_factor),min_dstretch)
    elif r > upper_bound:
        stretch_step_max = stretch_step_max*equil_factor
    return r

    
    


        

pwd = (!pwd)
pwd = pwd[0]

nwalkers = 500
active_box = nwalkers+1
nchains = 32
nbeads = 2
max_vol_per_atom = 15

ns_iterations = int(2e5)
walk_length = 50

aspect_ratio_limit = 0.8
angle_limit_deg = 55
angle_limit_rad = angle_limit_deg*np.pi/180


mdl.box_set_num_boxes(nwalkers+1) #nwalkers+2 if debugging
mdl.box_initialise()
mdl.box_set_pbc(1)
mdl.alkane_set_nchains(nchains) 
mdl.alkane_set_nbeads(nbeads)    
mdl.alkane_initialise()           
mdl.box_set_isotropic(1)
mdl.box_set_bypass_link_cells(1) # Bypass use of link cell algorithm for neighbour finding
mdl.box_set_use_verlet_list(0)   # Don't use Verlet lists either since CBMC moves quickly invalidate these

#random seeds
#np.random.seed(1)
#mdl.random_set_random_seed(1)

cell_matrix = 0.999*np.eye(3)*np.cbrt(nbeads*nchains*max_vol_per_atom)#*np.random.uniform(0,1)
for ibox in range(1,nwalkers+1):
    mdl.box_set_cell(ibox,cell_matrix)
    
    
#setting step sizes
mdl.alkane_set_dr_max(0.65)
mdl.alkane_set_dt_max(0.43)
mdl.alkane_set_dh_max(0.4)
mdl.alkane_set_dv_max(0.5)
shear_step_max = 0.5
stretch_step_max = 0.5

min_dr = 1e-10
min_dt = 1e-10
min_dh = 1e-10
min_dv = 1e-10
min_dshear = 1e-5
min_dstretch = 1e-5


#setting params for MC_run
move_types = ['box','translate', 'rotate', 'dihedral', 'shear', 'stretch']
ivol = 0; itrans = 1; irot = 2; idih = 3; ishear = 4; istr = 5
move_ratio = np.zeros(6)
move_ratio[ivol] = 10
move_ratio[itrans] = 3.0*nchains
move_ratio[irot] = 2.0*nchains
#move_ratio[idih] = 1.0*nchains
move_ratio[ishear] = 3
move_ratio[istr] = 3
#moves_prob = np.cumsum(moves_ratio)/np.sum(moves_ratio)
walk_length = 20


#constructing dictionaries which contain initial volumes
volumes = {}
start_volumes = []
for ibox in range(1,nwalkers+1):
    volumes[ibox], rate = MC_run(walk_length, move_ratio, ibox)
    start_volumes.append(volumes[ibox])
    
    
#overlap check
overlap_check = np.zeros(nwalkers)
for ibox in range(1,nwalkers+1):
    overlap_check[mdl.alkane_check_chain_overlap(ibox)]
# print(np.where(overlap_check!=0))


# main driver code


#removed_volumes = []


mc_adjust_interval = nwalkers//2 #for adjusting step sizes
low_acc_rate = 0.2
high_acc_rate = 0.5 

snapshots = 1000

vis_interval = ns_iterations//snapshots

filename = f"nested_sampling.{nwalkers}.1.extxyz"

try:
    os.remove(filename)
except:
    pass


energies_filename = f"nested_sampling.{nwalkers}.1.energies"

try:
    os.remove(energies_filename)
except:
    pass

energies_file = open(energies_filename, "w+")
energies_file.write(f"{nwalkers} {1} {5*nchains} {False} {nchains} \n")

#traj = ase.io.write("nestedsampling.extxyz", mode="a", fmt="extxyz")


# f = IntProgress(min=0, max=ns_iterations) 
# display(f) # display the bar

for i in range(ns_iterations):
    index_max = max(volumes, key=volumes.get)
    volume_limit = volumes[index_max]
    clone = np.random.randint(1,nwalkers+1)
    
    if i%vis_interval == 0:
        largest_config = mk_ase_config(index_max,nbeads,nchains)
        largest_config.wrap()
        ase.io.write(filename,largest_config, append = True)
#         traj.write(largest_config)


    
    if i%mc_adjust_interval == 0:
        rv = adjust_dv(clone,active_box,low_acc_rate,high_acc_rate, volume_limit)
        adjust_dr(clone,active_box,low_acc_rate,high_acc_rate)
        adjust_dt(clone,active_box,low_acc_rate,high_acc_rate)
#         adjust_dh(clone,active_box,low_acc_rate,high_acc_rate)
        adjust_dshear(clone,active_box,low_acc_rate,high_acc_rate)
        adjust_dstretch(clone,active_box,low_acc_rate,high_acc_rate)

    
    clone_walker(clone,active_box) #copies the ibox from first argument onto the second one.
    
    
    new_volume,_ = MC_run(walk_length, move_ratio, active_box,volume_limit)
    
    
    if new_volume-volume_limit<sys.float_info.epsilon:
        #removed_volumes.append(volumes[index_max])
        energies_file.write(f"{i} {volumes[index_max]} {volumes[index_max]} \n")
        volumes[index_max] = new_volume #replacing the highest volume walker
        clone_walker(active_box, index_max)
#     if i%100 == 0:
# #         print(i,volume_limit,rv)
#         pass

        
    
    

    f.value+=1

#traj.close()

# print("--- %s seconds ---" % (time.time() - start_time))


#overlap check
overlap_check = np.zeros(nwalkers)
for ibox in range(1,nwalkers+1):
    overlap_check[mdl.alkane_check_chain_overlap(ibox)]
# print(np.where(overlap_check!=0))
    