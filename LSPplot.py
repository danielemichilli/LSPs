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

from Parameters import *

def plot(idL,args):

  #-----------------------------------
  #Selection of values in the database
  #-----------------------------------

#  store = pd.HDFStore('SinlgePulses.hdf5','r')
#  print store
#  data = store[idL]
#  puls = store[idL+'_pulses']
#  store.close()
  
  
  puls = pd.read_hdf('SinlgePulses.hdf5',idL+'_pulses',where=['Pulse==0'])
  puls_rfi = pd.read_hdf('SinlgePulses.hdf5',idL+'_pulses',where=['(Pulse>0)&(Pulse<=args.rfi)'])
  
  
  
  #k = 4.1488078e3  #s
  #delay = np.float32(.5 * k * (F_MIN**-2 - F_MAX**-2))
  #puls.Time = puls.Time - puls.DM * delay
  #puls_rfi.Time = puls_rfi.Time + puls_rfi.DM * delay
  
  
  
  
  #puls=puls[puls.Sigma>6.5]

  #if not puls.shape[0]:
    #print 'The DataBase is empty.'
    #return
  
  puls = puls[puls.BEAM>12]
  puls = puls[(puls['DM']>args.dmlo) & (puls['DM']<args.dmhi) & (puls['Time']>args.tmlo) & (puls['Time']<args.tmhi) & (puls['SAP'].isin(args.sap)) & (puls['BEAM'].isin(args.beam))]
  
  puls_rfi = puls_rfi[puls_rfi.BEAM>12]
  puls_rfi = puls_rfi[(puls_rfi['DM']>args.dmlo) & (puls_rfi['DM']<args.dmhi) & (puls_rfi['Time']>args.tmlo) & (puls_rfi['Time']<args.tmhi) & (puls_rfi['SAP'].isin(args.sap)) & (puls_rfi['BEAM'].isin(args.beam))]
  
  data = pd.read_hdf('SinlgePulses.hdf5',idL,where=['Pulse=puls.index.tolist()'])
  
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
      #col = data.SAP *10 + (data.BEAM-13) *10./62.
      
      col = puls.SAP *10 + (puls.BEAM-13) *10./62.
      
      #col = puls.index.values.astype(float)
    
    else:
      col = data.Pulse.values
  
  else:
    col=u'r' 
    
    
  fig = plt.figure()

  
  ax1 = plt.subplot2grid((1,5),(0,0),colspan=3)
  ax2 = plt.subplot2grid((1,5),(0,3),sharey=ax1)
  ax3 = plt.subplot2grid((1,5),(0,4),sharey=ax1)

    
  ax1.scatter(data.Time, data.DM, facecolors='none', s=sig, c='k')

  ax1.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt=None, ecolor=col)
  
  ax1.scatter(puls_rfi.Time, puls_rfi.DM, s=20., c=u'k',marker='+')

  ax1.scatter(puls.Time, puls.DM, facecolors=col, s=150.)

  if args.c & (np.size(args.beam) + np.size(args.sap) > 2) & (puls.shape[0]):   #Testare che faccia tutto bene, sembra troppo robusto
    ticks = np.linspace(col.min(),col.max(),num=10)
    bar = plt.colorbar(ticks=ticks)
    bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
    bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
    bar.update_ticks
    bar.ax.xaxis.set_ticks_position('top')


  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([args.tmlo[0],args.tmhi[0],args.dmlo[0],args.dmhi[0]])
  
  ax2.hist(puls.DM.tolist(),bins=1000,orientation=u'horizontal')
  ax2.set_xlabel('Counts')

  ax3.scatter(puls.Sigma,puls.DM)
  ax3.set_xlabel('Sigma')
  
  fig.subplots_adjust(wspace=0)   
  for ax in [ax2, ax3]:
    plt.setp(ax.get_yticklabels(), visible=False)
    # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
    #ax.set_xticks(ax.get_xticks()[1:])  

  plt.show()
  


def RFI_channels(data):
  
  data.Time.hist(bins=3600)
  plt.show()

  plt.xlabel('Time (s)')
  plt.ylabel('Number of events')
  
  return