# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 20:57:46 2021

@author: brndn
"""

import numpy as np
import matplotlib.pyplot as plt

tau_min_param = 20
tau_max_param = 1000
shift = 0

v = np.linspace(-130,50,1000)

tau_s = tau_min_param + tau_max_param / ( np.exp((v+71.5-shift)/14.2) + np.exp(-(v+89-shift)/11.6) )

plt.figure()
plt.plot(v,tau_s)
plt.xlabel("Membrane Voltage (mV)")
plt.ylabel("Tau (ms)")

tau_max = np.amax(tau_s)
max_idx = np.argmax(tau_s)
tau_min = np.amin(tau_s)
min_idx = np.argmin(tau_s)
vmin = v[min_idx]
vmax = v[max_idx]
