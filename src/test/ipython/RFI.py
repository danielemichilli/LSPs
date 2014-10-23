import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import numpy as np

idL = os.path.basename(os.getcwd())
store = pd.HDFStore('SinlgePulses.hdf5','r')
data = store[idL]
puls = store[idL+'_pulses']
store.close()

puls = puls[puls.Pulse==0]
data = data.astype(np.float64)
gb = data[data.Pulse.isin(puls.index)].groupby('Pulse',sort=False)

puls1 = puls[(puls.DM<40)|(puls.DM>50)]

puls2 = puls[(puls.DM>43.46)&(puls.DM<43.48)]
#puls2 = puls[(puls.DM>=26.72)&(puls.DM<=26.82)]

data1 = data[data.Pulse.isin(puls1.index)]
data2 = data[data.Pulse.isin(puls2.index)]
gb1 = data1.groupby('Pulse',sort=False)
gb2 = data2.groupby('Pulse',sort=False)





float(len(puls1))/l1*100
float(len(puls2))/l2*100






float(len(b[b<]))/len(b)*100
float(len(a[a>]))/len(a)*100


a.hist(bins=100,alpha=.5,color='r')
b.hist(bins=100,alpha=.5,color='b')
plt.xlabel('Parameter value')
plt.ylabel('N')
plt.show()


plt.axis([,,0,])
plt.plot([,],[0,],color='g',linewidth=3.)