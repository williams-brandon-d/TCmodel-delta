# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 11:54:40 2021

@author: brndn
"""

from neuron import h 
import numpy as np
import pandas as pd
from neuron_sim_functions import *
    
ih_param_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
ih_param_filename = "Ih_model_values_CPZ.csv"

worksheet_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/"
worksheet_filename = "CPZ_TC0_Voltage_Simulations_rev1"

h.load_file("stdrun.hoc") # for run control
h.nrn_load_dll("C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/nrn_mech.dll")
h.load_file('C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/Fdelta.oc')

cell_names = ["TC[0]","TC[1]","RE[0]","RE[1]"]

norm_ghbar = 1e-5 # default max Ih conductance (S/cm^2)
h.tstop = 5000 # stop time in ms

trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1 = h.Vector(),h.Vector(),h.Vector(),h.Vector(),h.Vector()

trec.record(h._ref_t) # record time
vrec_tc0.record(h.TC[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_tc1.record(h.TC[1].soma[0](0.5)._ref_v) # record voltage from center
vrec_re0.record(h.RE[0].soma[0](0.5)._ref_v) # record voltage from center
vrec_re1.record(h.RE[1].soma[0](0.5)._ref_v) # record voltage from center

v_tc0_df = pd.DataFrame() # initialize dataframe for TC[0] voltage traces

ih_param_df = pd.read_csv(ih_param_path + ih_param_filename)

cd_control = getControlCD_CPZ(ih_param_df)
#cd_control_dict = getControlCDdict(ih_param_df)

print("cd_control = %.4g (pA/pF)" % cd_control)
print("")

for i in range(len(ih_param_df['Name'])):
    
    save_filename = ih_param_df['Name'].iloc[i]
    
    print("")
    print(save_filename)
    

   # set new ih params for both TC cells
    h.shift_iar =  75 + ih_param_df['v_half'].iloc[i]
    h.slope_iar = ih_param_df['slope_factor'].iloc[i]
    h.tau130_iar = ih_param_df['tau_fast'].iloc[i]
    
    # group_name = save_filename[:-2]
    # if  group_name.find('control') < 0:
    #     group_name = group_name + "control"   
    # print("group_name = %s" % group_name)
    # cd_control = cd_control_dict[group_name]
    # print("cd_control = %.4g (pA/pF)" % cd_control)
 
    
    for iCell in range(int(h.ncells)):    
        h.TC[iCell].kl.gmax = 0.004
        cd = ih_param_df['current_density'].iloc[i] 
        if iCell == 0: print("cd = %.4g (pA/pF)" % cd)
        h.TC[iCell].soma[0].ghbar_iar = norm_ghbar * cd / cd_control
        
    h.frecord_init()
    h.run()
    
    time,v_tc0,v_tc1,v_re0,v_re1 = vecToNumpy(trec,vrec_tc0,vrec_tc1,vrec_re0,vrec_re1)
    
    data_matrix = np.concatenate((time,v_tc0,v_tc1,v_re0,v_re1)).reshape(5,trec.size()).transpose()
    
    inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict = makeIBIdict(data_matrix,cell_names,printnum=1)
    
    plotSim(data_matrix,cell_names,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict,save_filename)
    
    saveFig(save_filename)
    
    saveData(save_filename,[data_matrix,inter_IBI_dict,inter_freq_dict,nAP_dict,mean_nAP_dict,time_AP1_dict,intra_IBI_dict,intra_freq_dict])
    
    if i == 0: v_tc0_df["time_ms"] = time
    
    v_tc0_df[save_filename] = v_tc0
    
#v_tc0_df.to_excel("%s%s.xlsx" % (worksheet_path,worksheet_filename))
