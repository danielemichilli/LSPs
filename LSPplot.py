#######################################################
#
# Written by Daniele Michilli
#
########################################################

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import numpy as np



def plot(idL,args):

  #-----------------------------------
  #Selection of values in the database
  #-----------------------------------

  store = pd.HDFStore('SinlgePulses.hdf5','r')
  data = store[idL]

  data=data[(data['DM']>args.dmlo) & (data['DM']<args.dmhi) & (data['Time']>args.tmlo) & (data['Time']<args.tmhi) & (data['SAP'].isin(args.sap)) & (data['BEAM'].isin(args.beam))]

  if args.s: 
    sig=(data.Sigma/5.)**3
  else:
    sig=8

  if args.l: 
    linewidths=data.Downfact*0.1
  else:
    linewidths=1
  
  if args.c: 
    #col = data.SAP *10 + (data.BEAM-13) *10./62.
  
    col = data.Pulse.astype(float)
  
  else:
    col='b'  

  print args.dmhi

  plt.scatter(data.Time, data.DM, facecolors='none', s=sig, c=col, cmap=mpl.cm.rainbow)

  if args.c:   #Testare che faccia tutto bene, sembra troppo robusto
    ticks = np.linspace(col.min(),col.max(),num=10)
    bar = plt.colorbar(ticks=ticks)
    bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
    bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
    bar.update_ticks
    bar.ax.xaxis.set_ticks_position('top')
  
  
  #confrontare plot e scatter: velocita e bellezza
  #plt.plot(data['Time'], data['DM'], 'ko', mfc='none', ms=2)

  plt.xlabel('Time (s)')
  plt.ylabel('DM (pc/cm3)')
  plt.axis([args.tmlo[0],args.tmhi[0],args.dmlo[0],args.dmhi[0]])

  plt.show()

  store.close()
