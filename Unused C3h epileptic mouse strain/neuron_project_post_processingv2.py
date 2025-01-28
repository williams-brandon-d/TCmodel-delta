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


def loadData(save_filename):
    data_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Data"
    data_filename = "%s.pckl" % save_filename   
    f = open(os.path.join(data_path,data_filename), 'rb')
    [data_matrix,_,inter_freq_dict,_,mean_nAP_dict,time_AP1_dict,_,intra_freq_dict] = pickle.load(f)
    f.close()
    return data_matrix,inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict

def plotBar(label_list,mean_list,xlabel_string,std_list):
    rev_label_list = label_list.copy()
    rev_label_list.reverse()
    rev_mean_list = mean_list.copy()
    rev_mean_list.reverse()
    rev_std_list = std_list.copy()
    rev_std_list.reverse()
    
    x_pos = [i for i, _ in enumerate(rev_label_list)]
    
    plt.figure()
    plt.barh(x_pos, rev_mean_list, color='blue', xerr=rev_std_list)
    plt.xlabel(xlabel_string,fontweight="bold")
    plt.yticks(x_pos, rev_label_list,fontweight="bold")
    plt.xticks(fontweight="bold")
    plt.show()

ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
ih_param_filename = "Ih_model_values_all.csv"
ih_param_df = pd.read_csv(os.path.join(ih_param_path,ih_param_filename))

key = "TC[0]"

group_names = ["bl6","c3hcontrol","c3hday1","c3hday7","c3hday25"]
cells_per_group = 7

inter_freq_mean = []
mean_nAP_mean = []
time_AP1_mean = []
intra_freq_mean = []
inter_freq_std = []
mean_nAP_std = []
time_AP1_std = []
intra_freq_std = []

for iGroup in range(len(group_names)):
    #fig, ax = plt.subplots(cells_per_group,1,figsize=(10,10))
    
    inter_freq_list = []
    mean_nAP_list = []
    time_AP1_list = []
    intra_freq_list = []
    
    for iCell in range(cells_per_group):
        
        save_filename = "%s_%d" % (group_names[iGroup],iCell+1)
        if save_filename == "c3hday1_3": continue # no spikes with parameters
        if save_filename == "c3hday7_7": continue
        
        data_matrix,inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict = loadData(save_filename)
        
        inter_freq_list.append(inter_freq_dict[key])
        mean_nAP_list.append(mean_nAP_dict[key])
        time_AP1_list.append(time_AP1_dict[key])
        intra_freq_list.append(intra_freq_dict[key])           
        
        # ax[iCell].plot(data_matrix[:,0],data_matrix[:,1])
        # #ax[iCell].set_title(" %s %s Interburst Freq = %.3g Hz, Time to 1st AP = %.3g ms, Intraburst Freq = %.3g Hz, APs/burst = %.3g" % (save_filename,key,inter_freq_dict[key],time_AP1_dict[key],intra_freq_dict[key],mean_nAP_dict[key]),fontweight="bold")
        # ax[iCell].set_title("%s" % save_filename,fontweight="bold")
        # ax[iCell].set_ylabel("Voltage (mV)",fontweight="bold")
        # ax[iCell].set_xlim(data_matrix[0,0],data_matrix[-1,0])
        # ax[iCell].set_ylim(-120,60)
        
        # #fontsize = 14
        # for tick in ax[iCell].xaxis.get_major_ticks():
        #     #tick.label1.set_fontsize(fontsize)
        #     tick.label1.set_fontweight('bold')
        # for tick in ax[iCell].yaxis.get_major_ticks():
        #     #tick.label1.set_fontsize(fontsize)
        #     tick.label1.set_fontweight('bold')
 
    inter_freq_mean.append(np.mean(inter_freq_list))
    mean_nAP_mean.append(np.mean(mean_nAP_list))
    time_AP1_mean.append(np.mean(time_AP1_list))
    intra_freq_mean.append(np.mean(intra_freq_list))
    
    inter_freq_std.append(np.std(inter_freq_list))
    mean_nAP_std.append(np.std(mean_nAP_list))
    time_AP1_std.append(np.std(time_AP1_list))
    intra_freq_std.append(np.std(intra_freq_list))
        
    # fig.tight_layout()
    # ax[-1].set_xlabel("Time (ms)",fontweight="bold")
    
plotBar(group_names,inter_freq_mean,"Interburst Frequency (Hz)",inter_freq_std)
plotBar(group_names,intra_freq_mean,"Intraburst Frequency (Hz)",intra_freq_std)
plotBar(group_names,time_AP1_mean,"Time to First AP (ms)",time_AP1_std)
plotBar(group_names,mean_nAP_mean,"APs/burst",mean_nAP_std)
