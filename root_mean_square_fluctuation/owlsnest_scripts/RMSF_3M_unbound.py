import sys, os
import mdtraj as md
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random
#----------------------------------------------------------------------------------
START = 200
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN12/curated_RUN12_p16925/RUN12CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r12_res = dict()
for i in range(0,700):
    M3_r12_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r12_res", M3_r12_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r12_res[j]
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
        np.save('./individual_means/M3_r12_res%d' %j, group1_avg)
#----------------------------------------------------------------------------------
START = 200
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN13/curated_RUN13_p16925/RUN13CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r13_res = dict()
for i in range(0,700):
    M3_r13_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r13_res", M3_r13_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r13_res[j]
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
        np.save('./individual_means/M3_r13_res%d' %j, group1_avg)
#----------------------------------------------------------------------------------
START = 200
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN14/curated_RUN14_p16925/RUN14CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r14_res = dict()
for i in range(0,700):
    M3_r14_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r14_res", M3_r14_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r14_res[j]
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
        np.save('./individual_means/M3_r14_res%d' %j, group1_avg)
#----------------------------------------------------------------------------------
START = 200
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN15/curated_RUN15_p16925/RUN15CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r15_res = dict()
for i in range(0,700):
    M3_r15_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r15_res", M3_r15_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r15_res[j]
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
        np.save('./individual_means/M3_r15_res%d' %j, group1_avg)
#----------------------------------------------------------------------------------
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN30/curated_RUN30_p16925/RUN30CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r30_res = dict()
for i in range(0,700):
    M3_r30_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r30_res", M3_r30_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r30_res[j]
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
        np.save('./individual_means/M3_r30_res%d' %j, group1_avg)

#----------------------------------------------------------------------------------
STRIDE = 210
STOP = 100000000
ref = md.load('../../gro_files/p16925-r31-c0.gro')         # ref
ref_ind = ref.topology.select("protein and name CA")             # residue 1
rmsf = dict()
for i in range(10000):
    try:
        t = md.load('../RUN31/curated_RUN31_p16925/RUN31CLONE%d.xtc'%i,top='../../gro_files/p16925-r31-c0.gro')
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
M3_r31_res = dict()
for i in range(0,700):
    M3_r31_res[i] = [item[i] for item in rmsf_list]
np.save("M3_r31_res", M3_r31_res)
rmsf_avg_grp_1 = []
for j in range(700):
    res_list = M3_r31_res[j]
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
        np.save('./individual_means/M3_r31_res%d' %j, group1_avg)

