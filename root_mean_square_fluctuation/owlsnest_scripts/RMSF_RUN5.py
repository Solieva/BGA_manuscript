import sys, os
import mdtraj as md
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random

START = 200
STRIDE = 210
STOP = 100000000
ref = md.load('../p14188-r6-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../curated_RUN5_p14188/RUN5CLONE%d.xtc'%i,top='../p14188-r6-c0.gro')
        t = t[START:STOP:STRIDE]              
        atom_indices1=t.topology.select("protein and name CA") 
        t.superpose(ref, atom_indices=atom_indices1, ref_atom_indices=ref_ind)
        avg_xyz = np.mean(t.xyz[:, atom_indices1, :], axis=0)
        rmsf[i] = np.sqrt(3*np.mean((t.xyz[:, atom_indices1, :] - avg_xyz)**2, axis=(0,2)))
    except (IOError):         
        pass
for i in range(10000): #removing all the arrays that had only 1 frame
    try: 
        if min(rmsf[i]) <= 0.0000000:
            del rmsf[i]
    except KeyError:
        pass 

new_rmsf = {k: v for k, v in rmsf.items() if pd.Series(v).notna().all()} 
rmsf_list = list(new_rmsf.values())
M2_r5_res = dict()
for i in range(0,700):
    M2_r5_res[i] = [item[i] for item in rmsf_list]


np.save("M2_r5_res", M2_r5_res)


rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M2_r5_res[j]
    rmsf_avg_rmsf = { i : res_list[i] for i in range(0, len(res_list) ) }
    rmsf_avg_grp_1.append(rmsf_avg_rmsf)


group1_avg = []
for j in range(700):
    group1_avg = []
    for i in range(1000):  # boot strapping
        d = rmsf_avg_grp_1[j]
        keys = random.sample(list(d), 100)
        rmsf_values = [d[k] for k in keys]
        rmsf_list = rmsf_values
        n = len(rmsf_list)                        
        def average(a,n):                  
            sum = 0                            
            for i in range(n):                 
                sum += a[i]                       
            return sum/n;                      
        average_rmsf = average(rmsf_list, n)
        avg1 = np.average(average_rmsf)
        group1_avg.append(avg1)
        group1_avg.sort(reverse=False)
        np.save('./individual_means/M2_r5_res%d' %j, group1_avg)


