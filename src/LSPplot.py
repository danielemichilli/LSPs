#############################
#
# LOTAAS Single Pulse plots
#
# Written by Daniele Michilli
#
#############################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tarfile
import os
import matplotlib as mpl

from Parameters import *

def plot(puls,puls_rfi,meta_data,top_candidates,best_pulses,color=True,store=False,data=pd.DataFrame()):
  
  #puls_rfi = puls[puls.Sigma<=6.]
  #puls = puls[puls.Sigma>6.]
  
  col = puls.Sigma
  if color: 
    cmap = plt.get_cmap('gist_heat_r')
    fill = u'b'
    square = u'g'
  else: 
    cmap = plt.get_cmap('Greys')
    fill = u'k'
    square = u'k'
    
  sig = (top_candidates.Sigma/1.5)**3

  fig = plt.figure()
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  ax5 = plt.subplot2grid((3,4),(0,3))

  ax1.scatter(puls_rfi.Time, puls_rfi.DM, s=5., c=u'k',marker='+',linewidths=[0.4,])
  
 
  if not puls.empty:
    if not data.empty: 
      ax1.scatter(data.Time, data.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])
      
    ax1.scatter(puls.Time, puls.DM, c=col, s=20., cmap=cmap,linewidths=[0.,],vmin=5,vmax=10)
    ax1.plot([0,3600],[40.5,40.5],'k--')
    ax1.plot([0,3600],[141.7,141.7],'k--')
    
    if not top_candidates.empty: ax1.scatter(top_candidates.Time,top_candidates.DM,s=sig,linewidths=[0.,],c=fill,marker='*')
    if not best_pulses.empty: ax1.scatter(best_pulses.Time,best_pulses.DM,s=sig,linewidths=[1.,],marker='s',facecolors='none',edgecolor=square)
    
    ax1.set_yscale('log')
        
    mpl.rc('font', size=5)
    for i in range(0,top_candidates.shape[0]):
      ax1.annotate(i,xy=(top_candidates.Time.iloc[i],top_candidates.DM.iloc[i]*1.15),horizontalalignment='center',verticalalignment='bottom')
    
    if len(puls.DM.unique())>1:
      ax2.hist(puls.Sigma.tolist(),bins=100,histtype='step',color='k')
      ax2.set_xlabel('SNR')
      ax2.set_ylabel('Counts')
      ax2.set_yscale('log')
      
      hist = ax3.hist(puls.DM.tolist(),bins=300,histtype='stepfilled',color=u'k')
      ax3.set_xscale('log')
      ax3.set_xlabel('DM (pc/cm3)')
      ax3.set_ylabel('Counts')
      ax3.set_xlim(5,550)
      ax3.plot([40.5,40.5],[0,hist[0].max()],'k--')
      ax3.plot([141.7,141.7],[0,hist[0].max()],'k--')
      
    ax4.scatter(puls.DM,puls.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=5,vmax=10)
    ax4.scatter(top_candidates.DM,top_candidates.Sigma,s=15.,linewidths=[0.,],c=fill,marker='*')
    ax4.scatter(best_pulses.DM,best_pulses.Sigma,s=15.,linewidths=[1.,],c='b',marker=u's',facecolors='none',edgecolor=square)
    ax4.set_xscale('log')
    ax4.set_ylabel('SNR')
    ax4.set_xlabel('DM (pc/cm3)')
    ax4.axis([5,550,puls.Sigma.min(),puls.Sigma.max()+3.])
    ax4.plot([40.5,40.5],[0,best_pulses.Sigma.max()],'k--')
    ax4.plot([141.7,141.7],[0,best_pulses.Sigma.max()],'k--')
      
    mpl.rc('font', size=3.5)
    for i in range(0,top_candidates.shape[0]):
      ax4.annotate(i,xy=(top_candidates.DM.iloc[i]*1.15,top_candidates.Sigma.iloc[i]),horizontalalignment='left',verticalalignment='center')
    
    mpl.rc('font', size=5)
    ax5.axis([0,10,0,7])
    ax5.annotate('File: '+meta_data.File.iloc[0], xy=(0,6))
    ax5.annotate('Telescope: '+meta_data.Telescope.iloc[0], xy=(0,5))
    ax5.annotate('Instrument: '+meta_data.Instrument.iloc[0], xy=(0,4))
    ax5.annotate('RA: '+meta_data.RA.iloc[0], xy=(0,3))
    ax5.annotate('DEC: '+meta_data.DEC.iloc[0], xy=(0,2))
    ax5.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,1))
    ax5.axis('off')

  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,5,550])
  
  ax1.tick_params(which='both',direction='out')
  ax2.tick_params(which='both',direction='out')
  ax3.tick_params(which='both',direction='out')
  ax4.tick_params(which='both',direction='out')
  ax5.tick_params(which='both',direction='out')
  
  fig.tight_layout()

  if store:
    mpl.rc('font',size=5)
    plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  else: 
    plt.show()
  
  plt.clf()
  return
  


def sp(top_candidates,data,meta_data,size=True,store=False):
  
  fig = plt.figure()
  
  for i in range(0,top_candidates.shape[0]):
  
    puls = top_candidates.iloc[i]
    events = data[data.Pulse==puls.name]
    
    if size: sig=(events.Sigma/5.)**4
    else: sig=8
  
    ax = plt.subplot2grid((2,5),(i/5,i%5))
    ax.scatter(events.Time, events.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])  
    ax.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt='none', ecolor='r')
    #ax.set_title = "Pulse "+str(i)+" (DM "+str(puls.DM)+")"
    #ax.set_xlabel('Time (s)')
    #ax.set_ylabel('DM (pc/cm3)')
    
    fig.tight_layout()
    
  #fig.text(0.5, 0.04, 'Time (s)', ha='center', va='center')
  #fig.text(0.06, 0.5, 'DM (pc/cm3)', ha='center', va='center', rotation='vertical')
    
  if store:
    if not top_candidates.empty:
      mpl.rc('font',size=5)
      plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  else: plt.show()
  
  plt.clf()
  return


def obs_top_candidates(top_candidates,best_pulses,color=True,size=True,store=False,incoherent=False): #top_candidates di tutti i beams
  
  if color:
    if incoherent:
      col_top = top_candidates.SAP
      if not best_pulses.empty: col_best = best_pulses.SAP
    else:
      col_top = top_candidates.SAP *10 + (top_candidates.BEAM-13) /61. * 10.
      if not best_pulses.empty: col_best = best_pulses.SAP *10 + (best_pulses.BEAM-13) /61. * 10.
  else: 
    col_top = u'r' 
    col_best = u'r'
  
  if size: 
    sig_top = (top_candidates.Sigma/6.)**4
    if not best_pulses.empty: sig_best = (top_candidates.Sigma/6.)**4
  else: sig=100.
    
  fig = plt.figure()
  plt.title("Best Candidates")
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=3,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))

  main_plt = ax1.scatter(top_candidates.Time,top_candidates.DM,s=sig_top,linewidths=[0.,],c=col_top)
  
  dim = len(top_candidates.SAP.unique())+len(top_candidates.BEAM.unique())-1
  
  if color & (dim>1):
    if incoherent:
      ticks = np.linspace(col_top.min(),col_top.max(),num=3)
      bar = plt.colorbar(mappable=main_plt,ticks=ticks,ax=ax1)
      bar.set_ticklabels(['{0:.0f}'.format(int(t)) for t in ticks])
      bar.ax.set_xlabel('sap',ha='left',labelpad=10)
      bar.update_ticks
      bar.ax.xaxis.set_ticks_position('top')      
    else:
      ticks = np.linspace(col_top.min(),col_top.max(),num=10)
      bar = plt.colorbar(mappable=main_plt,ticks=ticks,ax=ax1)
      bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10*6.+13) for t in ticks])
      bar.ax.set_xlabel('sap, beam',ha='left',labelpad=10)
      bar.update_ticks
      bar.ax.xaxis.set_ticks_position('top')
    
    
  if not best_pulses.empty: ax1.scatter(best_pulses.Time,best_pulses.DM,s=sig_best,linewidths=[1.,],marker='s',facecolors='none',c=col_best)
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,5,550])
  ax1.set_yscale('log')

  if not top_candidates.empty:
    if len(top_candidates.DM.unique())>1:
      ax3.hist(top_candidates.DM.tolist(),bins=300,histtype='stepfilled',color=u'k')
      ax3.set_xscale('log')
      ax3.set_xlabel('DM (pc/cm3)')
      ax3.set_ylabel('Counts')
      ax3.set_xlim(5,550)
    
      ax2.hist(top_candidates.Sigma.tolist(),bins=100,histtype='step',color='k')
      ax2.set_xlabel('SNR')
      ax2.set_ylabel('Counts')
      ax2.set_yscale('log')
    
    ax4.scatter(top_candidates.DM,top_candidates.Sigma,c=u'k',s=3.,linewidths=[0.,],vmin=5,vmax=10)
    ax4.scatter(top_candidates.DM,top_candidates.Sigma,s=15.,linewidths=[0.,],c=u'k',marker='*')
    ax4.scatter(best_pulses.DM,best_pulses.Sigma,s=15.,linewidths=[1.,],c=u'b',marker=u's',facecolors='none',edgecolor=u'g')
    ax4.set_xscale('log')
    ax4.set_ylabel('SNR')
    ax4.set_xlabel('DM (pc/cm3)')
    ax4.axis([5,550,top_candidates.Sigma.min(),top_candidates.Sigma.max()+3.])
  
  ax1.tick_params(which='both',direction='out')
  ax2.tick_params(which='both',direction='out')
  ax3.tick_params(which='both',direction='out')
  ax4.tick_params(which='both',direction='out')
  
  fig.tight_layout()
  
  if store:
    mpl.rc('font',size=5)
    plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  else: 
    plt.show()
  
  plt.clf()
  return
    
      
