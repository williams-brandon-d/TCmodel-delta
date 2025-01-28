# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 11:54:40 2021

@author: brndn
"""
#c3hday1_3, c3hday7_7  : no spikes

from neuron import h 
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

def plotSim(data_matrix,titles,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict):
    nPlots = data_matrix.shape[1] - 1
    fig, ax = plt.subplots(nPlots,1,figsize=(10,10))

    for iAxis in range(len(ax)):
        key = titles[iAxis]
        ax[iAxis].plot(data_matrix[:,0],data_matrix[:,iAxis+1])
        #ax[iAxis].set_title("%s Interburst Freq = %.3g Hz, Time to 1st AP = %.3g ms, Intraburst Freq = %.3g Hz, APs/burst = %.3g" % (key,inter_freq_dict[key],time_AP1_dict[key],intra_freq_dict[key],mean_nAP_dict[key]))
        ax[iAxis].set_ylabel("Voltage (mV)")
        ax[iAxis].set_xlim(data_matrix[0,0],data_matrix[-1,0])
        ax[iAxis].set_ylim(-120,60)
        
    fig.tight_layout()
    ax[-1].set_xlabel("Time (ms)")

def calculateIBIs(data,time,plotnum=0):
    max_peaks,_ = find_peaks(data,height=0)
    min_peaks,_ = find_peaks(-1*data,height=-1*np.mean(data))
    # print("max_peaks = ",max_peaks)
    time_AP1 = time[max_peaks[0]]
    temp_max_peaks = max_peaks
    burst_AP1 = []
    burst_nAPs = []
    intraburst_IBIs = []
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
        inter_freq_dict[key] = 1000/np.mean(inter_IBI_dict[key])
        intra_freq_dict[key] = 1000/np.mean(intra_IBI_dict[key])
        mean_nAP_dict[key] = np.mean(nAP_dict[key])
        if printnum and iCell == 0: print("%s interburst_freq = %.3g Hz" % (key,inter_freq_dict[key]))
    return inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict

    
ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
ih_param_filename = "Ih_model_values_all.csv"
ih_param_df = pd.read_csv(os.path.join(ih_param_path,ih_param_filename))

h.load_file("stdrun.hoc") # for run control

h.nrn_load_dll("C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/nrn_mech.dll")

h.load_file('C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/Fdelta.oc')

trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1 = h.Vector(),h.Vector(),h.Vector(),h.Vector(),h.Vector()

trec.record(h._ref_t) # record time
vrec_tc0.record(h.TC[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_tc1.record(h.TC[1].soma[0](0.5)._ref_v) # record voltage from center
vrec_re0.record(h.RE[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_re1.record(h.RE[1].soma[0](0.5)._ref_v) # record voltage from center

cell_names = ["TC[0]","TC[1]","RE[0]","RE[1]"]

norm_ghbar = 1e-5 # default max Ih conductance (S/cm^2)

h.tstop = 2000 # stop time in ms

for i in range(len(ih_param_df['Name'])):
    
    save_filename = ih_param_df['Name'].iloc[i]
    if save_filename == "c3hday1_3": continue # no spikes with parameters
    if save_filename == "c3hday7_7": continue

   # set new ih params for both TC cells
    h.shift_iar =  75 + ih_param_df['v_half'].iloc[i]
    h.slope_iar = ih_param_df['slope_factor'].iloc[i]
    #h.taumax_iar = ih_param_df['tau_slow'].iloc[i]
    h.tau130_iar = ih_param_df['tau_fast'].iloc[i]
 
    
    for iCell in range(int(h.ncells)):    
        h.TC[iCell].kl.gmax = 0.004
        h.TC[iCell].soma[0].ghbar_iar = norm_ghbar*ih_param_df['current_density'].iloc[i]/ih_param_df['current_density'].iloc[0]

    h.frecord_init()
    h.run()
    
    time,v_tc0,v_tc1,v_re0,v_re1 = vecToNumpy(trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1)
    
    data_matrix = np.concatenate((time,v_tc0,v_tc1,v_re0,v_re1)).reshape(5,trec.size()).transpose()
            
    printnum = 1
    if printnum: print(ih_param_df['Name'].iloc[i])
    
    inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict = makeIBIdict(data_matrix,cell_names,printnum)
    
    plotSim(data_matrix,cell_names,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict)
    
    saveFig(save_filename)
    
    saveData(save_filename,[data_matrix,inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict])

