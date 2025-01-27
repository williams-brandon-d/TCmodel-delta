# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 14:37:29 2021

@author: brndn
"""

from neuron import h 
import numpy as np
from neuron_sim_functions import makeIBIdict,getControlCD,saveData

cd_range = [5, 20] # (pA/pF)
vhalf_range = [-110, -50] # (mV)
slope_range = [5, 15] # (mV^-1)
tau_range = [100, 5000] # (ms)

nPoints = 9

vectorDict = {}
vectorDict['cd'] = np.linspace(cd_range[0],cd_range[1],16)
vectorDict['vhalf'] = np.linspace(vhalf_range[0],vhalf_range[1],16)
vectorDict['slope'] = np.linspace(slope_range[0],slope_range[1],nPoints)
vectorDict['tau'] = np.linspace(tau_range[0],tau_range[1],nPoints)

h.load_file("stdrun.hoc") # for run control
h.nrn_load_dll("C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/nrn_mech.dll")
h.load_file('C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/Fdelta.oc')

cell_names = ["TC[0]"]

norm_ghbar = 1e-5 # default max Ih conductance (S/cm^2)
h.tstop = 5000 # stop time in ms

trec,vrec_tc0 = h.Vector(),h.Vector()

trec.record(h._ref_t) # record time
vrec_tc0.record(h.TC[0].soma[0](0.5)._ref_v) # record voltage from center

cd_control = getControlCD()

nCDPoints = len(vectorDict['cd'])
nVhalfPoints = len(vectorDict['vhalf'])
nSlopePoints = len(vectorDict['slope'])
nTauPoints = len(vectorDict['tau'])

tc0_inter_freq_array = np.zeros((nCDPoints,nVhalfPoints,nSlopePoints,nTauPoints),dtype = np.float64())

print("cd_control = %.4g (pA/pF)" % cd_control)
print("")

for iCD in range(nCDPoints):
    cd = vectorDict['cd'][iCD]
    for iVhalf in range(nVhalfPoints):
        vhalf = vectorDict['vhalf'][iVhalf]
        for iSlope in range(nSlopePoints):
            slope = vectorDict['slope'][iSlope]
            for iTau in range(nTauPoints):
                tau = vectorDict['tau'][iTau]
                
               # set new ih params for both TC cells
                h.shift_iar =  75 + vhalf
                h.slope_iar = slope
                h.tau130_iar = tau
                for iCell in range(int(h.ncells)):    
                    h.TC[iCell].kl.gmax = 0.004
                    #if iCell == 0: print("cd = %.4g (pA/pF)" % cd)
                    h.TC[iCell].soma[0].ghbar_iar = norm_ghbar * cd / cd_control
                
                print("")
                string = "cd = %.3g, vhalf = %.3g, slope = %.3g, tau = %.3g" % (cd,vhalf,slope,tau)
                print(string)
                
                h.frecord_init()
                h.run()
                
                time = trec.as_numpy()
                v_tc0 = vrec_tc0.as_numpy()
                
                data_matrix = np.concatenate((time,v_tc0)).reshape(2,trec.size()).transpose()
                
                _,inter_freq_dict,_,_,_,_,_ = makeIBIdict(data_matrix,cell_names,printnum=0)
                
                if np.array(inter_freq_dict["TC[0]"]).size > 0:
                    tc0_inter_freq_array[iCD,iVhalf,iSlope,iTau] = np.mean(inter_freq_dict["TC[0]"])
                
                print("burst freq = %.3g Hz" % tc0_inter_freq_array[iCD,iVhalf,iSlope,iTau])
                
                #plotSim(data_matrix,cell_names,inter_freq_dict,time_AP1_dict,intra_freq_dict,mean_nAP_dict,string)
                
saveData("parameter_sweep_nPoints%d" % nPoints,[tc0_inter_freq_array,vectorDict])