# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 16:50:42 2021

@author: brndn
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle

def loadParameterSweepData(save_filename):
    data_path = "C:/Users/brndn/OneDrive/Desktop/White Lab/NEURON Project/Data/"
    data_filename = "%s.pckl" % save_filename   
    f = open(data_path+data_filename, 'rb')
    [inter_freq_array,vectorDict] = pickle.load(f)
    f.close()
    return inter_freq_array,vectorDict

def minmax(data):
    return np.amin(data),np.amax(data)

matplotlib.rcParams['axes.unicode_minus'] = True

nPoints = 9

save_filename = "parameter_sweep_nPoints%d" % nPoints
inter_freq_array,vectorDict = loadParameterSweepData(save_filename)

cdMin,cdMax = minmax(vectorDict['cd'])
vhalfMin,vhalfMax = minmax(vectorDict['vhalf'])
slopeMin,slopeMax = minmax(vectorDict['slope'])
tauMin,tauMax = minmax(vectorDict['tau'])

dummy_tau = np.linspace(tauMin,tauMax,5)
plot_tau_vector = np.concatenate((dummy_tau[0],vectorDict['tau'][1],dummy_tau[1:]),axis=None)

plot_slope_vector = np.linspace(slopeMin,slopeMax,5)

tau_indices = np.where(np.in1d(vectorDict['tau'], plot_tau_vector))[0]
slope_indices = np.where(np.in1d(vectorDict['slope'], plot_slope_vector))[0]

vhalf_ticks = [-100,-50]
vhalf_ticklabels = ["−100  ","−50"]
cd_ticks = [10, 20]

max_freq = np.amax(inter_freq_array)
im_extent = [ vhalfMin, vhalfMax, cdMin, cdMax ]
aspect_ratio = ( vhalfMax - vhalfMin ) / ( cdMax - cdMin )
nRows = len(plot_slope_vector)
nCols = len(plot_tau_vector)

fig, axes = plt.subplots(nrows=nRows,ncols=nCols)

for iSlope in range(nRows):
    slope_index = slope_indices[iSlope]
    for iTau in range(nCols):
        tau_index = tau_indices[iTau]
        ax = axes[-(iSlope+1),iTau]
        im = ax.imshow(inter_freq_array[:,:,slope_index,tau_index], origin ='lower', extent =im_extent,
                   aspect=aspect_ratio, cmap='hot', interpolation='nearest', vmin=0, vmax=max_freq)
        
        if iSlope == nRows - 1:
            ax.set_title("%.3d" % plot_tau_vector[iTau],fontsize=10,fontweight='bold')
            
        if iTau == nCols - 1:
            ax.text(1.05, 0.5, "%.3g" % plot_slope_vector[iSlope],
                rotation=0, size=10, weight='bold',
                ha='left', va='center', transform=ax.transAxes)
            
        if iSlope == 0 and iTau == 0:
            ax.set_xlabel("$\mathbf{V_{0.5}}$ (mV)",fontsize=10,fontweight='bold')
            ax.set_ylabel("Current Density\n (pA/pF)",fontsize=8,fontweight='bold')
            ax.set_xticks(vhalf_ticks)
            ax.set_yticks(cd_ticks)
            ax.set_xticklabels(vhalf_ticklabels, fontsize=9,fontweight='bold')
            ax.set_yticklabels(cd_ticks, fontsize=9,fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        else:
            ax.set_axis_off()
        
cbar = fig.colorbar(im, ax=axes.ravel().tolist(),pad=0.25)
cbar.ax.set_title("Interburst\n Frequency (Hz)",fontweight='bold')

#fig.suptitle("Tau fast (ms)",fontweight='bold',ha='right')
#fig.suptitle(r'$\tau_{fast}$ (ms)',fontweight='bold',ha='right')
fig.text(0.3,0.95,r'$\tau_{fast}$ (ms)', size=12, weight='bold',ha='left', va='center')

#fig.text(0.68,0.5,"Slope Factor $\mathbf{(mV^{-1}})$",rotation=90, size=12, weight='bold',ha='left', va='center')
fig.text(0.68,0.5,"Slope Factor $\mathbf{(mV)}$",rotation=90, size=12, weight='bold',ha='left', va='center')

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.show()

#fig.savefig('%s_final.tif' % save_filename,dpi=300)
fig.savefig('Oniani22_fig14_final.tif',dpi=300)