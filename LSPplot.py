#############################
#
# LOTAAS Single Pulse plots
#
# Written by Daniele Michilli
#
#############################

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
  puls = store[idL+'_pulses']
  store.close()
  
  #data = data[(data['DM']>args.dmlo) & (data['DM']<args.dmhi) & (data['Time']>args.tmlo) & (data['Time']<args.tmhi) & (data['SAP'].isin(args.sap)) & (data['BEAM'].isin(args.beam))]
  
  #puls = puls[puls.BEAM<16]
  
  #puls=puls[puls.SAP==0]
  
  #data = data[data.BEAM<16]
  
  #data=data[data.SAP==0]
  
  #data.sort('BEAM')
  
  #data.groupby(['SAP','BEAM','Pulse'],sort=False).DM.count()>3
  
  puls = puls[puls.Pulse>0]
  
  puls = puls[puls.BEAM>12]
  
  puls = puls[(puls['DM']>args.dmlo) & (puls['DM']<args.dmhi) & (puls['Time']>args.tmlo) & (puls['Time']<args.tmhi)]
  
  if not puls.shape[0]:
    print 'The DataBase is empty.'
    return
  
  data = data[data.Pulse.isin(puls.index)]

  if args.s: 
    sig=(data.Sigma/5.)**3
  else:
    sig=8

  if args.l: 
    linewidths=data.Downfact*0.1
  else:
    linewidths=1
    
  if args.c:
    
    if np.size(args.beam) + np.size(args.sap) > 2:
      col = data.SAP *10 + (data.BEAM-13) *10./62.
      
      #col = puls.SAP *10 + (puls.BEAM-13) *10./62.
      
      #col = puls.index.values.astype(float)
    
    else:
      col = data.Pulse.values
  
  else:
    col='b'  
  
  plt.scatter(data.Time, data.DM, facecolors='none', s=sig, c=col)#, cmap=mpl.cm.rainbow)
  
  if args.c & (np.size(args.beam) + np.size(args.sap) > 2):   #Testare che faccia tutto bene, sembra troppo robusto
    ticks = np.linspace(col.min(),col.max(),num=10)
    bar = plt.colorbar(ticks=ticks)
    bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
    bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
    bar.update_ticks
    bar.ax.xaxis.set_ticks_position('top')
    
    
  plt.scatter(puls.Time, puls.DM, facecolors='none', s=50.)#, c=col)
  
  plt.errorbar(puls.Time, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt=None, c=col) #ecolor=puls.index.values.astype(float), 
    
  
  #plt.scatter(puls.Time, puls.DM, marker='+')


  
  #confrontare plot e scatter: velocita e bellezza
  #plt.plot(data['Time'], data['DM'], 'ko', mfc='none', ms=2)

  plt.xlabel('Time (s)')
  plt.ylabel('DM (pc/cm3)')
  plt.axis([args.tmlo[0],args.tmhi[0],args.dmlo[0],args.dmhi[0]])

  plt.show()
  


def RFI_channels(data):
  
  data.Time.hist(bins=3600)
  plt.show()

  plt.xlabel('Time (s)')
  plt.ylabel('Number of events')
  
  return