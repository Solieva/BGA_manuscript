# import packages:
import sys, os
import mdtraj as md
import numpy as np
print("packages loaded")
# define the calculation:
class SASA:
    def __init__(self,name):
        self.name = name

    def sasa(structure_file,clones,path):
        ref = md.load(structure_file) # structure_file = structure file to superpose trajectories to (ex: '.gro file')
        START = 200
        STRIDE = 210
        STOP = 100000000

        for c in range(clones): # for each trajectory 
            print(c)
            try:
                traj = md.load(path+'%d.xtc'%c, top=structure_file)  # load trajectory
                traj = traj[START:STOP:STRIDE]
                atom_indices1=traj.topology.select("protein and name CA") # trajectory indices for superimposing
                ref_ind = ref.topology.select("protein and name CA") # structure file indices for superimposing
                traj.superpose(ref, atom_indices=atom_indices1, ref_atom_indices=ref_ind) # superpose each trajectory  
                atom_indices2=traj.topology.select("protein") # all atoms in protein
                small_traj = traj.atom_slice(atom_indices=atom_indices2) # slice
                res = md.shrake_rupley(small_traj, probe_radius=0.14, n_sphere_points=960, mode='atom') # sasa calculation
                sasa1.append(res)
            except (IOError):
                pass

clones = 494
run_number_1M = []
run_number_2M = [5,6,7]
run_number_3M = []
run_number_4M = []
run_numbers = [run_number_1M,run_number_2M,run_number_3M,run_number_4M]

structure_path = '../../../curated_p14188_trajectories/'
structure_numbers = [0,6,0,0]
molarity_number = ["1","2","3","4"]

for i in [1]:
    structure_file = structure_path+'p14188-r%d-c0.gro'%structure_numbers[i]
    print(structure_file)
    for j in range(len(run_numbers[i])):
        path = "../../../curated_p14188_trajectories/curated_RUN%d_p14188/RUN%dCLONE"%(run_numbers[i][j],run_numbers[i][j])
        sasa1 = []
        SASA.sasa(structure_file,clones,path)
        sasa_per_traj = []
        np.save("SASA_data_files/SASA_%sM_run%d_unbound"%(molarity_number[i],run_numbers[i][j]),sasa1) 
