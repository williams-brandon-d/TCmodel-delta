# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 13:28:29 2021

@author: brndn
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.signal import find_peaks
import pickle

def vecToNumpy(trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1):
    time = trec.as_numpy()
    v_tc0 = vrec_tc0.as_numpy()
    v_tc1 = vrec_tc1.as_numpy()
    v_re0 = vrec_re0.as_numpy()
    v_re1 = vrec_re1.as_numpy()
    return time,v_tc0,v_tc1,v_re0,v_re1

def saveData(save_filename,obj):
    data_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Data"
    data_filename = "%s.pckl" % save_filename    
    f = open(os.path.join(data_path,data_filename), 'wb')
    pickle.dump(obj, f)
    f.close()
    
def saveFig(save_filename):
    fig_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Figures"
    fig_filename = "%s.png" % save_filename   
    plt.savefig(os.path.join(fig_path,fig_filename))

def plotSim(data_matrix,titles,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict,plot_title):
    nPlots = data_matrix.shape[1] - 1
    fig, ax = plt.subplots(nPlots,1,figsize=(10,10))

    for iAxis in range(len(ax)):
        #key = titles[iAxis]
        ax[iAxis].plot(data_matrix[:,0],data_matrix[:,iAxis+1])
        #ax[iAxis].set_title("%s Interburst Freq = %.3g Hz, Time to 1st AP = %.3g ms, Intraburst Freq = %.3g Hz, APs/burst = %.3g" % (key,inter_freq_dict[key],time_AP1_dict[key],intra_freq_dict[key],mean_nAP_dict[key]))
        ax[iAxis].set_ylabel("Voltage (mV)")
        ax[iAxis].set_xlim(data_matrix[0,0],data_matrix[-1,0])
        ax[iAxis].set_ylim(-120,60)
        
    fig.tight_layout()
    ax[-1].set_xlabel("Time (ms)")
    fig.suptitle(plot_title)

def calculateIBIs(data,time,plotnum=0,printnum=0):
    interburst_IBIs =[]
    burst_nAPs = []
    time_AP1 = []
    intraburst_IBIs = []
    
    max_peaks,_ = find_peaks(data,height=0)
    
    if len(max_peaks) == 0:
        if printnum: print("No Spikes Detected.")
        return interburst_IBIs, burst_nAPs, time_AP1, intraburst_IBIs
    
    min_peaks,_ = find_peaks(-1*data,height=-1*np.mean(data))
    # print("max_peaks = ",max_peaks)
    time_AP1 = time[max_peaks[0]]
    temp_max_peaks = max_peaks
    burst_AP1 = []

    for iii in range(len(max_peaks)):
        if (temp_max_peaks.size == 0): break
        first_max_peak = temp_max_peaks[0]
        burst_AP1.append(first_max_peak)   
        temp_min_peaks = min_peaks[min_peaks > first_max_peak]
        if (temp_min_peaks.size == 0): break
        next_min_peak = temp_min_peaks[0]
        burst_APs = temp_max_peaks[temp_max_peaks < next_min_peak]
        if len(burst_APs) > 1:            
            intraburst_IBIs.append(time[burst_APs[1]] - time[burst_APs[0]])
            burst_nAPs.append(len(burst_APs))
        temp_max_peaks = max_peaks[max_peaks > next_min_peak]
        # print("first_max_peak = %g" % first_max_peak)
        # print("next_min_peak = %g" % next_min_peak)
        # print("temp_max_peaks =",temp_max_peaks)    
    if (plotnum == 1):
        plt.figure()
        plt.plot(data)
        plt.scatter(burst_AP1,data[burst_AP1],marker="o",color="r")
        # plt.scatter(data_matrix[max_peaks,0],data_matrix[max_peaks,col],marker="o",color="r")
        # plt.scatter(data_matrix[min_peaks,0],data_matrix[min_peaks,col],marker='o',color='k')
    interburst_IBIs = np.diff(time[burst_AP1])
    return interburst_IBIs, burst_nAPs, time_AP1, intraburst_IBIs

def makeIBIdict(data_matrix,dict_key,printnum=1):
    inter_IBI_dict = {}
    intra_IBI_dict = {}
    nAP_dict = {}
    mean_nAP_dict = {}
    inter_freq_dict = {}
    intra_freq_dict = {}
    time_AP1_dict = {}
    for iCell in range(data_matrix.shape[1]-1):
        key = dict_key[iCell]
        inter_IBI_dict[key],nAP_dict[key],time_AP1_dict[key],intra_IBI_dict[key] = calculateIBIs(data_matrix[:,iCell+1],data_matrix[:,0])
        if len(inter_IBI_dict[key]) == 0:
            inter_freq_dict[key] = []
            if printnum and iCell == 0: print("%s Only One Burst Detected." % key)

        else:
            inter_freq_dict[key] = 1000/np.mean(inter_IBI_dict[key])
            if printnum and iCell == 0: print("%s interburst_freq = %.3g Hz" % (key,inter_freq_dict[key]))

        if len(intra_IBI_dict[key]) == 0:
            intra_freq_dict[key] = []
        else:
            intra_freq_dict[key] = 1000/np.mean(intra_IBI_dict[key])
            
        if len(nAP_dict[key]) == 0:
            mean_nAP_dict[key] = []
        else:
            mean_nAP_dict[key] = np.mean(nAP_dict[key])
    return inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict

# control for IFNalpha and ILbeta simulations
def getControlCD():
    ih_control_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
    ih_control_filename = "Ih_model_values.csv"
    fullfile = '%s%s' % (ih_control_path,ih_control_filename)
    ih_param_df = pd.read_csv(fullfile)
    cd_control = ih_param_df['current_density'].iloc[0] # average bl6 current density
    return cd_control

# # control for CPZ simulations   
def getControlCD_CPZ(ih_param_df): # calculate average current density of all controls
    indices = [i for i, elem in enumerate(ih_param_df['Name']) if "control" in elem]
    return np.mean(ih_param_df['current_density'][indices])

# not used as any controls    
# def getControlCDdict(ih_param_df):
#     cd_control_dict = {}
#     group_name_list = []
    
#     #find unique group names
#     for iCell in range(len(ih_param_df['Name'])):
#         save_filename = ih_param_df['Name'].iloc[iCell]
#         group_name = save_filename[:-2]
#         if  group_name.find('control') > 0:
#             group_name = group_name[:-len('control')]
#         group_name_list.append(group_name)
        
#     unique_group_names = np.unique(np.array(group_name_list))
    
#     # calculate average current density for each control group
#     for iGroup in range(len(unique_group_names)):
#         unique_group_name = unique_group_names[iGroup]
#         control_group_name = unique_group_name + "control"
#         indices = [i for i, elem in enumerate(ih_param_df['Name']) if control_group_name in elem]
#         cd_list = ih_param_df['current_density']
#         control_cd_list = cd_list[indices]
#         cd_control_dict[control_group_name] = np.mean(control_cd_list)
        
#         # print("control_group_name = %s" % control_group_name)
#         # print("control_cd_list")
#         # print(control_cd_list)
#         # print("average_control = %.3g (pA/pF)" % cd_control_dict[group_name])

#     return cd_control_dict
