# import packages:
import sys, os
import mdtraj as md
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random

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
run_number_1M = [5,6,7,26,27]
run_number_2M = [38,39]
run_number_3M = [13,14,15,30,31]
run_number_4M = [21,22,23,34,35]
run_numbers = [run_number_1M,run_number_2M,run_number_3M,run_number_4M]

structure_path = '../../gro_files/'
structure_numbers = [27,39,31,35]
molarity_number = ["1","2","3","4"]

for i in range(len(run_numbers)):
    structure_file = structure_path+'p16925-r%d-c0.gro'%structure_numbers[i]
    #print(structure_file)
    for j in range(len(run_numbers[i])):
        path = "../RUN%d/curated_RUN%d_p16925/RUN%dCLONE"%(run_numbers[i][j],run_numbers[i][j],run_numbers[i][j])
        sasa1 = []
        SASA.sasa(structure_file,clones,path)
        sasa_per_traj = []
        np.save("SASA_data_files/SASA_%sM_run%d_unbound"%(molarity_number[i],run_numbers[i][j]),sasa1) 
