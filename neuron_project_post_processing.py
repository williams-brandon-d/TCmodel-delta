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

def plotBar(label_list,data_list,xlabel_string):
    rev_label_list = label_list.copy()
    rev_label_list.reverse()
    rev_data_list = data_list.copy()
    rev_data_list.reverse()
    x_pos = [i for i, _ in enumerate(rev_label_list)]
    
    plt.figure()
    plt.barh(x_pos, rev_data_list, color='blue')
    plt.xlabel(xlabel_string,fontweight="bold")
    plt.yticks(x_pos, rev_label_list,fontweight="bold")
    plt.xticks(fontweight="bold")
    plt.show()

ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
ih_param_filename = "Ih_model_values_all.csv"
ih_param_df = pd.read_csv(os.path.join(ih_param_path,ih_param_filename))

key = "TC[0]"

nPlots = len(ih_param_df['Name'])

fig, ax = plt.subplots(nPlots,1,figsize=(10,10))

inter_freq_list = []
mean_nAP_list = []
time_AP1_list = []
intra_freq_list = []
names_list = []

for i in range(nPlots):
    
    save_filename = ih_param_df['Name'].iloc[i]
    data_matrix,inter_freq_dict,mean_nAP_dict,time_AP1_dict,intra_freq_dict = loadData(save_filename)
    
    ax[i].plot(data_matrix[:,0],data_matrix[:,1])
    #ax[i].set_title(" %s %s Interburst Freq = %.3g Hz, Time to 1st AP = %.3g ms, Intraburst Freq = %.3g Hz, APs/burst = %.3g" % (save_filename,key,inter_freq_dict[key],time_AP1_dict[key],intra_freq_dict[key],mean_nAP_dict[key]),fontweight="bold")
    ax[i].set_title("%s" % save_filename,fontweight="bold")
    ax[i].set_ylabel("Voltage (mV)",fontweight="bold")
    ax[i].set_xlim(data_matrix[0,0],data_matrix[-1,0])
    ax[i].set_ylim(-120,60)
    
    #fontsize = 14
    for tick in ax[i].xaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax[i].yaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')

    inter_freq_list.append(inter_freq_dict[key])
    mean_nAP_list.append(mean_nAP_dict[key])
    time_AP1_list.append(time_AP1_dict[key])
    intra_freq_list.append(intra_freq_dict[key])
    names_list.append(save_filename)
    
fig.tight_layout()
ax[-1].set_xlabel("Time (ms)",fontweight="bold")

plotBar(names_list,inter_freq_list,"Interburst Frequency (Hz)")
plotBar(names_list,intra_freq_list,"Intraburst Frequency (Hz)")
plotBar(names_list,time_AP1_list,"Time to First AP (ms)")
plotBar(names_list,mean_nAP_list,"APs/burst")
