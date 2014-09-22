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
import tarfile

from Parameters import *

def plot(idL,puls,puls_rfi,data,meta_data,top_candidates,size=False,color=False,store=False):

  
  if size: 
    sig=(data.Sigma/5.)**4
  else:
    sig=8

  if color:
    
    if data.SAP.unique().shape[0] + data.BEAM.unique().shape[0] > 2:
      #col = data.SAP *10 + (data.BEAM-13) *10./62.
      
      col = puls.SAP *10 + (puls.BEAM-13) *10./62.
      
      #col = puls.index.values.astype(float)
    
    else:
      col = puls.index.values

  else:
    col=u'r' 
    
    
  fig = plt.figure()

  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,1))
  ax3 = plt.subplot2grid((3,4),(0,2))
  ax4 = plt.subplot2grid((3,4),(0,0))
  ax5 = plt.subplot2grid((3,4),(0,3))

    
  ax1.scatter(puls_rfi.Time, puls_rfi.DM, s=20., c=u'k',marker='+')
  
  if not puls.empty: 

    ax1.scatter(data.Time, data.DM, facecolors='none', s=sig, c='k')

    ax1.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt=None, ecolor='k')#ecolor=col)
  
    ax1.scatter(puls.Time, puls.DM, facecolors=col, s=150.)

    if color & (puls.shape[0]):   #Testare che faccia tutto bene, sembra troppo robusto
      ticks = np.linspace(col.min(),col.max(),num=10)
      bar = plt.colorbar(ticks=ticks)
      bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
      bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
      bar.update_ticks
      bar.ax.xaxis.set_ticks_position('top')

    ax2.hist(puls.DM.tolist(),bins=1000,histtype='stepfilled',color='k')
    ax2.set_xlabel('DM (pc/cm3)')
    ax2.set_ylabel('Counts')

    ax3.scatter(puls.Sigma,puls.DM,c='k')
    ax3.set_xlabel('SNR')
    ax3.set_ylabel('DM (pc/cm3)')
    ax3.axis([puls.Sigma.min(),puls.Sigma.max()+1.,0,550])
    
    ax4.hist(puls.Sigma.tolist(),bins=1000,histtype='stepfilled',color='k')
    ax4.set_xlabel('SNR')
    ax4.set_ylabel('Counts')
    
    print meta_data
    
    ax5.axis([0,10,0,6])
    ax5.annotate('File: '+str(meta_data[0]), xy=(1,5))
    ax5.annotate('Telescope: '+meta_data[1], xy=(1,4))
    ax5.annotate('RA: '+meta_data[2], xy=(1,3))
    ax5.annotate('DEC: '+meta_data[3], xy=(1,2))
    ax5.annotate('Epoch (MJD): '+meta_data[4], xy=(1,1))
    ax5.axis('off')

  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,0,550])
  
  #fig.subplots_adjust(wspace=0)   
  #for ax in [ax2, ax3]:
    #plt.setp(ax.get_yticklabels(), visible=False)

  if store:
    store = 1    
    
  else: plt.show()
  
  return
  


def store(idL,puls,puls_rfi,data,size=False,color=False):

  
  if size: 
    sig=(data.Sigma/5.)**4
  else:
    sig=8

  if color:
    
    if data.SAP.unique().shape[0] + data.BEAM.unique().shape[0] > 2:
      #col = data.SAP *10 + (data.BEAM-13) *10./62.
      
      col = puls.SAP *10 + (puls.BEAM-13) *10./62.
      
      #col = puls.index.values.astype(float)
    
    else:
      col = puls.index.values

  else:
    col=u'r' 
    
    
  fig = plt.figure()
  
  ax1 = plt.subplot2grid((2,4),(1,0),colspan=3)
  ax2 = plt.subplot2grid((2,4),(0,1))
  ax3 = plt.subplot2grid((2,4),(0,2))
  ax4 = plt.subplot2grid((2,4),(0,0))
    
  ax1.scatter(puls_rfi.Time, puls_rfi.DM, s=20., c=u'k',marker='+')
  
  if not puls.empty: 

    ax1.scatter(data.Time, data.DM, facecolors='none', s=sig, c='k')

    ax1.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt=None, ecolor='k')#ecolor=col)
  
    ax1.scatter(puls.Time, puls.DM, facecolors=col, s=150.)

    if color & (puls.shape[0]):   #Testare che faccia tutto bene, sembra troppo robusto
      ticks = np.linspace(col.min(),col.max(),num=10)
      bar = plt.colorbar(ticks=ticks)
      bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
      bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
      bar.update_ticks
      bar.ax.xaxis.set_ticks_position('top')

    ax2.hist(puls.DM.tolist(),bins=1000,histtype='stepfilled',color='k')
    ax2.set_xlabel('DM (pc/cm3)')
    ax2.set_ylabel('Counts')

    ax3.scatter(puls.Sigma,puls.DM,c='k')
    ax3.set_xlabel('SNR')
    ax3.set_ylabel('DM (pc/cm3)')
    ax3.axis([puls.Sigma.min(),puls.Sigma.max()+1.,0,550])
    
    ax4.hist(puls.Sigma.tolist(),bins=1000,histtype='stepfilled',color='k')
    ax4.set_xlabel('SNR')
    ax4.set_ylabel('Counts')

  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,0,550])
  
  
