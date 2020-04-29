"""
sublemental material, Fig. S1a

TODO: rename files without date and nbar->rhobar
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

foldername = "/project/theorie/h/Henrik.Weyer/2019-2020-th-Masterarbeit/Coarsening-letter/peak-plateau-crossover-many-peaks/"

#parameters
dl = 0.5
Dm = 1
Dc = 10000
#parameters of the single runs: [rhobar, number of run, domain length]
parray = [[[1.5,1,200000],[1.5,2,200000],[1.5,3,200000],[1.5,4,200000]],
          [[3,1,200000],[3,2,200000]],
           [[6,1,200000],[6,2,200000]]]
#scaling of time in the theory solution scalef=2/(rhobar*mu*nu)
scalef = [2.3,1.4,0.7]
mu = 1/0.75

#numerical solution of the coarsening rate from etaStat(M)_continuation.nb, coarsening master curve
ts = np.array([28.6696, 34.4526, 41.2306, 49.1751, 58.4873, 69.4029, 82.1986, 97.1995, 114.788, 135.413, 159.605, 187.989, 221.302, 260.414, 306.356, 360.349, 423.838, 498.543, 586.507, 690.164, 812.419, 956.742, 1127.29, 1329.06, 1568.06, 1851.51, 2188.17, 2588.64, 3065.81, 3635.39, 4316.59, 5132.98, 6113.61, 7294.37, 8719.84, 10445.6, 12541.2, 15094.3, 18215.8, 22046.5, 26766.9, 32609., 39873.7, 48953., 60362.2, 74783., 93124.8, 116611., 146903., 186275., 237874., 306100., 397168., 519964., 687343., 918145., 1.24034*10**6, 1.69607*10**6, 2.34976*10**6, 3.30142*10**6, 4.70897*10**6, 6.82608*10**6, 1.00677*10**7, 1.51259*10**7, 2.31787*10**7, 3.62742*10**7, 5.80549*10**7, 9.5155*10**7, 1.59964*10**8, 2.76239*10**8, 4.90828*10**8, 8.98867*10**8, 1.69967*10**9, 3.32469*10**9, 6.74097*10**9])
masses = 200*np.array([1.32407, 1.40252, 1.48563, 1.57366, 1.6669, 1.76567, 1.87029, 1.98112, 2.09851, 2.22285, 2.35456, 2.49408, 2.64186, 2.7984, 2.96422, 3.13986, 3.32591, 3.52298, 3.73173, 3.95285, 4.18707, 4.43517, 4.69797, 4.97634, 5.27121, 5.58354, 5.91439, 6.26484, 6.63606, 7.02927, 7.44578, 7.88697, 8.3543, 8.84932, 9.37368, 9.9291, 10.5174, 11.1406, 11.8008, 12.5, 13.2407, 14.0252, 14.8563, 15.7366, 16.669, 17.6567, 18.7029, 19.8112, 20.9851, 22.2285, 23.5456, 24.9408, 26.4186, 27.984, 29.6422, 31.3986, 33.2591, 35.2298, 37.3173, 39.5285, 41.8707, 44.3517, 46.9797, 49.7634, 52.7121, 55.8354, 59.1439, 62.6484, 66.3606, 70.2927, 74.4578, 78.8697, 83.543, 88.4932, 93.7368])


#plot
f, ax1 = plt.subplots(1, 1, figsize=(10, 12))

#sweep over different rhobars
for k,ps in enumerate(parray):
    #collect the trajectories of the peak numbers for the different runs
    peakns=[]
    #collect the times
    tts=[]
    #added-up length of the domains for one rhobar
    Ltot = 0
    
    rhobar = ps[0][0]

    #sweep over single runs
    for i,params in enumerate(ps):
        seed = params[1]
        #domain length
        L = params[2]
        Ltot += L
        
        #import data
        ef = h5py.File(foldername+"20200220-2cMcRD-rec-det-L-{:d}-dx-{:.1f}-Dc-10000-nbar-{:.1f}-seed-{:d}-pbc".format(L,dl,rhobar,seed).replace(".","_")+".h5","r")
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
        ax1.plot(et[amin:],dl*L/peakn[amin:], linewidth=1, color="grey")   
        tts.append(et[amin:])
        peakns.append(peakn[amin:])
        ef.close()

    #calculate mean trajectory of all runs for fixed rhobar
    
    peaknsn=[]
    #align time of trajectories
    if len(tts)>1:
        lenl = np.zeros(len(tts))
        lenu = np.zeros(len(tts))
        for i,li in enumerate(tts):
            lenl[i] = np.min(li)
            lenu[i] = np.max(li)
        shl = np.max(lenl)
        shu = np.min(lenu)
        tt = tts[0][np.where(tts[0]==shl)[0][0]:np.where(tts[0]==shu)[0][0]]
        for i,li in enumerate(tts):
            peaknsn.append(peakns[i][np.where(tts[i]==shl)[0][0]:np.where(tts[i]==shu)[0][0]])
    else:
        tt = tts[0]
    #plot mean
    ax1.plot(tt,dl*Ltot/np.sum(peaknsn, axis=0), linewidth=3, color="k")   
    
    #theory
    
    #plot power-law
    amax = np.argwhere(tt>3e4)[0][0]
    norm = np.argwhere(tt==1e3)[0][0]
    if k == 0:
        ax1.plot(tt[:amax],dl*L/peaknsn[-1][norm]*1.1*(tt[:amax]/1e3)**(3/8), linewidth=2, label = "~ t^(3/8)")
    #plot numerical coarsening rate
    ax1.plot(ts*scalef[k],1/(mu*rhobar)*masses, linewidth=2, linestyle="--", label="rhobar = {:.1f}".format(rhobar))

#plot
ax1.set_title('recruitment - detachment kinetics, 1D', fontsize=26)
ax1.set_ylabel('average peak separation', fontsize=26)
ax1.set_xlabel('time', fontsize=26)
ax1.set_xlim(1e1,1e9)
ax1.set_ylim(2e1,1e4)
ax1.legend(fontsize=20, loc="upper left")
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.setp(ax1.get_yticklabels(), fontsize=20)
plt.setp(ax1.get_xticklabels(), fontsize=20)

plt.show(f)