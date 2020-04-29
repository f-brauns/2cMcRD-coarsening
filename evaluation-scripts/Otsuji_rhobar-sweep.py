# -*- coding: utf-8 -*-
"""
sublemental material, Fig. S2

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
#from layout import colorsad as colors
#from layout import linetypesad as linetypes
from scipy.signal import find_peaks

foldername = ""

dl = 0.1

#parameters of single runs: [rhobar, number of run, domain length]
parray = [[[3,1,20000],[3,2,20000]],
          [[10,1,20000],[10,2,20000],[10,3,20000]]]

#plot
f, ax1 = plt.subplots(1, 1, figsize=(14, 10))

#sweep over different rhobars
for k,ps in enumerate(parray):
    #length of the domains added up for one rhobar
    Ltot = 0
    
    rhobar = ps[0][0]
    
    #sweep over single runs
    for i,params in enumerate(ps):
        seed = params[1]
        #domain length
        L = params[2]
        Ltot += L
        
        #import data
        ef = h5py.File(foldername+"2cMcRD-Otsuji-L-{:d}-dx-{:.1f}-rhobar-{:.1f}-seed-{:d}-pbc".format(L,dl,rhobar,seed).replace(".","_")+".h5","r")
        ed = ef["data"]
        et = ef["times"]
        times = et[:]

        #sweep over times and count the peaks
        peakn = np.zeros(len(et))
        for j,t in enumerate(et):
            peaks,_ = find_peaks(ed[j], prominence = 1)
            peakn[j] = len(peaks)
        
        #starting point where peaks can be discriminated
        amin = np.argmax(peakn)
        #plot single run
        ax1.plot(et[amin:],dl*L/peakn[amin:], linewidth=1, color="k")   
        ef.close()

#theory
ax1.plot([100,3e7],[2,2*(3e7/100)**(1/3)], linewidth = 2)
    
ax1.set_title("Otsuji's model 2, 1D", fontsize=26)
ax1.set_ylabel('average peak separation', fontsize=26)
ax1.set_xlabel('time', fontsize=26)
ax1.set_xlim(10,1e10)
ax1.set_ylim(5e-1,2e2)
ax1.legend(fontsize=10, loc=0)
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.setp(ax1.get_yticklabels(), fontsize=20)
plt.setp(ax1.get_xticklabels(), fontsize=20)

plt.show(f)