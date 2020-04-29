# -*- coding: utf-8 -*-
"""
Convert Comsol txt export to hdf5

We use here m,c instead of u,v as names for the densities, i.e. Dm,Dc instead of Du,Dv.
We use n instead of rho as name for the total density in the meta data.

Comsol export settings:
Data export of density rho
File type: Text
Points to evaluate in: Regular grid
Data format: Grid
Number of X points: choose appropriately
Include header checked
"""
import numpy as np
import h5py

#parameters, adopt appropriately
dim = 1
L = 200000
dx = [0.5]
xs = np.arange(dx[0]/2,L,dx[0])
gridsize = [len(xs)]

rhobar = 1.5
seed = 1
Dm = 1
Dc = 10000
kon = 1
kr = 1
Kd = 1

#metadata in hdf5 file
meta = {"PDE": "recruitment-detachment kinetics, peak-mesa crossover due to sloped FBS",
        "domain dimension": dim,
        "structure": "[t,x]",
        "dim": [0,len(xs)],
        "ddim": [0,dx[0]],
        "tmax": 0,
        "L": L,
        "nbar": rhobar,
        "Dm": Dm,
        "Dc": Dc,
        "kon": kon,
        "kr": kr,
        "Kd": Kd,
        "m_initial": "rhobar*(1+0.01*rand(x))",
        "c_initial": "0",
        "f(m,c)": "(kon+kr*m)*c-m/(Kd+m)"}

grid = xs

times = []
#open txt file, i.e. the comsol output, adopt appropriately
f = open("2cMcRD-rec-det-Dc-{:d}-rhobar-{:.1f}-L-{:d}-pbc-seed-{:d}".format(Dc,rhobar,L,seed).replace(".","_")+".txt")

#find the time points of the data
for i,line in enumerate(f):
    if "@" in line:
        linesplit = line.split()
        times = np.append(times, float(line.split(" ")[3].replace("t=","").replace(",","").replace("E","e")))
   
tmax = times[-1]
dt1 = times[1]-times[0]
dt2 = times[2]-times[1]
if dt1 == dt2:
    dt = dt1
else:
    dt = 0
nt = len(times)
meta["tmax"] = tmax
meta["ddim"][0] = dt
meta["dim"][0] = nt

#load data
sol = np.zeros((nt,len(xs)),dtype=float)
f.seek(0,0)
for line in f:
    if "@" in line:
        print("@ line", line)
        break

for it in range(nt):
    ll = f.readline()
    sol[it,:]=np.array([float(w) for w in ll.split()])                   
    f.readline()
    f.readline()
    
f.close()

#create export file
ef = h5py.File("2cMcRD-rec-det-L-{:d}-dx-{:.1f}-Dc-{:d}-rhobar-{:.1f}-seed-{:d}-pbc".format(L,dx[0],Dc,rhobar,seed).replace(".","_")+".h5","w-")

#export
ef.create_dataset("times", data=times, compression="gzip")
ef.create_dataset("grid", data=grid, compression="gzip")
ed = ef.create_dataset("data", data=sol, compression="gzip")
#add attributes
for ke, va in meta.items():
    ed.attrs[ke] = va
ef.close()