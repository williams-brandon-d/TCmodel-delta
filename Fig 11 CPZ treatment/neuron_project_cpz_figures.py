# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:46:19 2021
Post-processing of NEURON Simulations

@author: brndn
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from scipy import stats

def loadData(save_filename):
    data_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Data/"
    data_filename = "%s.pckl" % save_filename   
    f = open(os.path.join(data_path,data_filename), 'rb')
    [data_matrix,_,inter_freq_dict,_,mean_nAP_dict,time_AP1_dict,_,intra_freq_dict] = pickle.load(f)
    f.close()
    return data_matrix,inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict

def plotBar(label_list,mean_list,xlabel_string,sem_list):
    rev_label_list = label_list.copy()
    rev_label_list.reverse()
    rev_mean_list = mean_list.copy()
    rev_mean_list.reverse()
    rev_sem_list = sem_list.copy()
    rev_sem_list.reverse()
    
    x_pos = [i for i, _ in enumerate(rev_label_list)]
    
    plt.figure()
    plt.barh(x_pos, rev_mean_list, color='blue', xerr=rev_sem_list)
    plt.xlabel(xlabel_string,fontweight="bold")
    plt.yticks(x_pos, rev_label_list,fontweight="bold")
    plt.xticks(fontweight="bold")
    plt.show()

def findUniqueGroups(ih_param_df):
    group_name_list = [] # find unique group names
    for iCell in range(len(ih_param_df['Name'])):
        cell_name = ih_param_df['Name'].iloc[iCell]
        group_name = cell_name[:-2] # remove cell number from name
        if  group_name.find('control') > 0:
            group_name = group_name[:-len('control')] # remove control from name
        group_name_list.append(group_name) 
    unique_names = np.unique(np.array(group_name_list)).tolist()
    unique_names.sort() # sorts normally by alphabetical order
    unique_names.sort(key=len) # sorts by descending length
    return unique_names

ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
ih_param_filename = "Ih_model_values_CPZ.csv"
#ih_param_filename = "Ih_model_values_IFNa_ILb.csv"

group_names = ["CPZday1control","CPZday1","CPZday7control","CPZday7","CPZday25control","CPZday25"]
#group_names = ["BL6_IFNalpha","BL6_IL1beta"]

data_key = "TC[0]"

ih_param_df = pd.read_csv(ih_param_path+ih_param_filename)

cpzDict = {}

inter_freq_group = {}
mean_nAP_group = {}
time_AP1_group = {}
intra_freq_group = {}

inter_freq_mean = []
mean_nAP_mean = []
time_AP1_mean = []
intra_freq_mean = []

inter_freq_sem = []
mean_nAP_sem = []
time_AP1_sem = []
intra_freq_sem = []


for iGroup in range(len(group_names)):
    group_name = group_names[iGroup]  
    
    inter_freq_group[group_name] = []
    mean_nAP_group[group_name] = []
    time_AP1_group[group_name] = []
    intra_freq_group[group_name] = []
    
    # find number of cells in group
    list_of_names = ih_param_df['Name']
    nCells = 0
    for iName in range(len(list_of_names)):
        name = ih_param_df['Name'].iloc[iName]
        compare = int(name[:-2] == group_names[iGroup])
        nCells = nCells + compare

    for iCell in range(nCells):
        
        save_filename = "%s_%d" % (group_names[iGroup],iCell+1)
        
        data_matrix,inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict = loadData(save_filename)
        
        inter_freq_group[group_name].append(inter_freq_dict[data_key])
        mean_nAP_group[group_name].append(mean_nAP_dict[data_key])
        time_AP1_group[group_name].append(time_AP1_dict[data_key])
        intra_freq_group[group_name].append(intra_freq_dict[data_key])           
 
    inter_freq_mean.append(np.mean(inter_freq_group[group_name]))
    mean_nAP_mean.append(np.mean(mean_nAP_group[group_name]))
    time_AP1_mean.append(np.mean(time_AP1_group[group_name]))
    intra_freq_mean.append(np.mean(intra_freq_group[group_name]))
    
    inter_freq_sem.append(stats.sem(inter_freq_group[group_name]))
    mean_nAP_sem.append(stats.sem(mean_nAP_group[group_name]))
    time_AP1_sem.append(stats.sem(time_AP1_group[group_name]))
    intra_freq_sem.append(stats.sem(intra_freq_group[group_name]))
    
#unique_names = findUniqueGroups(ih_param_df)
  
plotBar(group_names,inter_freq_mean,"Interburst Frequency (Hz)",inter_freq_sem)
plotBar(group_names,intra_freq_mean,"Intraburst Frequency (Hz)",intra_freq_sem)
plotBar(group_names,time_AP1_mean,"Time to First AP (ms)",time_AP1_sem)
plotBar(group_names,mean_nAP_mean,"APs/burst",mean_nAP_sem)
