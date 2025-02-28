# Oniani T, Vinnenberg L, Chaudhary R, Schreiber JA, Riske K, Williams B, Pape H-C, White JA, Junker A, Seebohm G, et al. Effects of Axonal Demyelination, Inflammatory Cytokines and Divalent Cation Chelators on Thalamic HCN Channels and Oscillatory Bursting. International Journal of Molecular Sciences. 2022; 23(11):6285. https://doi.org/10.3390/ijms23116285

## NEURON Simulation Methods

Simulations were conducted within the NEURON simulation environment (Hines and Carnevale 2001) by modifying a thalamocortical (TC) 
and thalamic reticular network model (Destexhe et al 1996). The original model can be downloaded at http://senselab.med.yale.edu/ModelDB/ 
and accessed using model number 3343. The model simulates delta oscillations by connecting 2 spontaneously pacemaking TC cells and 
2 thalamic reticular neurons (NRT) via AMPA and GABA<sub>A+B</sub> synapses. Four I<sub>h</sub> parameters were modified based on experimental results: 
the maximum I<sub>h</sub> conductance (_g<sub>hbar</sub>_), the voltage of half maximal activation (_V<sub>0.5</sub>_), the slope of the voltage activation curve (_k_), 
and the time constant at -130 mV (_τ<sub>fast</sub>_). The percent change in maximum I<sub>h</sub> conductance was found by dividing the current density 
measured in each cell by the average of the combined control groups. The maximum I<sub>h</sub> conductance was normalized to the default conductance 
in the model (10 μS/cm<sup>2</sup>). The voltage of half maximal activation was shifted from the default of the model (-75 mV) to match the measured 
value in each cell. The shift in voltage activation curve was equally applied to the voltage-dependent time constant curve. The voltage dependent 
time constant curve was then scaled to match the measured time constant at -130 mV. The I<sub>h</sub> parameters were equally modified in both TC cells,
while parameters in NRT cells were unaltered. The maximum potassium leak conductance was set to 0.004 μS for both TC cells. All simulations were 
5 seconds long and the temperature was set to 36 °C. 
