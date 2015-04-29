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

def plot(gb):
  pulses = gb[0]
  rfi = gb[1]
  meta_data = gb[2]
  events = gb[3]
  folder = gb[4]
  obs = os.path.basename(folder)
  sap = gb[0].SAP.iloc[0]
  beam = gb[0].BEAM.iloc[0]
  
  plt.clf()

  if beam == 12:
    sp_shape(pulses.head(10),events,'{}/sp/SAP{}_BEAM{}/top_candidates(0-9).png'.format(folder,sap,beam),obs)
    sp_shape(pulses.iloc[10:20],events,'{}/sp/SAP{}_BEAM{}/top_candidates(10-19).png'.format(folder,sap,beam),obs)
    sp_shape(pulses.iloc[20:30],events,'{}/sp/SAP{}_BEAM{}/top_candidates(20-29).png'.format(folder,sap,beam),obs)
    plt.clf()
    sp_plot(pulses.iloc[30:],rfi,meta_data,pulses.head(30),sap,beam,'{}/sp/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
    
  else:
    sp_shape(pulses.head(10),events,'{}/sp/SAP{}_BEAM{}/top_candidates.png'.format(folder,sap,beam),obs)
    plt.clf()
    sp_plot(pulses.iloc[10:],rfi,meta_data,pulses.head(10),sap,beam,'{}/sp/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
  
  return



def sp_plot(pulses,rfi,meta_data,top_candidates,sap,beam,store,events=pd.DataFrame()):

  col = pulses.Sigma
  cmap = plt.get_cmap('gist_heat_r')
  fill = u'b'
  square = u'g'
    
  sig_top = (top_candidates.Sigma/1.5)**3

  fig = plt.figure()
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  ax5 = plt.subplot2grid((3,4),(0,3))

  ax1.scatter(pulses.Time, pulses.DM, c=col, s=20., cmap=cmap, linewidths=[0.,], vmin=5, vmax=10)    
  ax1.scatter(rfi.Time, rfi.DM, s=5., c=u'k', marker='+', linewidths=[0.4,])

  if not events.empty: 
    ax1.scatter(events.Time, events.DM, facecolors='none', s=sig, c='k', linewidths=[0.5,])
      
  ax1.plot([0,3600],[40.48,40.48],'k--')
  ax1.plot([0,3600],[141.68,141.68],'k--')
    
  if not top_candidates.empty: ax1.scatter(top_candidates.Time, top_candidates.DM, s=sig_top, linewidths=[0.,], c=fill, marker='*')

  ax1.set_yscale('log')
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,5,550])
  
  mpl.rc('font', size=5)
  for i in range(0,top_candidates.shape[0]):
    if top_candidates.Time.iloc[i]>100:
      ax1.annotate(i,xy=(top_candidates.Time.iloc[i]-60,top_candidates.DM.iloc[i]),horizontalalignment='right',verticalalignment='center')

  if len(pulses.DM.unique())>1:
    ax2.hist(pulses.Sigma.tolist(),bins=100,histtype='step',color='k')
    ax2.set_xlabel('SNR')
    ax2.set_ylabel('Counts')
    ax2.set_yscale('log')
  
    hist = ax3.hist(pulses.DM.tolist(),bins=300,histtype='stepfilled',color=u'k')
    ax3.set_xscale('log')
    ax3.set_xlabel('DM (pc/cm3)')
    ax3.set_ylabel('Counts')
    ax3.set_xlim(5,550)
    ax3.plot([40.5,40.5],[0,hist[0].max()+10],'k--')
    ax3.plot([141.68,141.68],[0,hist[0].max()+10],'k--')
  
  ax4.scatter(pulses.DM,pulses.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=5,vmax=10)
  ax4.scatter(top_candidates.DM,top_candidates.Sigma,s=15.,linewidths=[0.,],c=fill,marker='*')
  ax4.set_xscale('log')
  ax4.set_ylabel('SNR')
  ax4.set_xlabel('DM (pc/cm3)')
  limit = top_candidates.Sigma.max()
  ax4.axis([5,550,pulses.Sigma.min(),limit+3.])
  ax4.plot([40.5,40.5],[0,limit+3.],'k--')
  ax4.plot([141.68,141.68],[0,limit+3.],'k--')
  mpl.rc('font', size=3.5)
  for i in range(0,top_candidates.shape[0]):
    ax4.annotate(i,xy=(top_candidates.DM.iloc[i]/1.15,top_candidates.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center')

  mpl.rc('font', size=5)
  ax5.axis([0,10,0,7])
  ax5.annotate('File: '+meta_data.File.iloc[0], xy=(0,6))
  ax5.annotate('Telescope: '+meta_data.Telescope.iloc[0], xy=(0,5))
  ax5.annotate('Instrument: '+meta_data.Instrument.iloc[0], xy=(0,4))
  ax5.annotate('RA: '+meta_data.RA.iloc[0], xy=(0,3))
  ax5.annotate('DEC: '+meta_data.DEC.iloc[0], xy=(0,2))
  ax5.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,1))
  ax5.axis('off')

  ax1.tick_params(which='both',direction='out')
  ax2.tick_params(which='both',direction='out')
  ax3.tick_params(which='both',direction='out')
  ax4.tick_params(which='both',direction='out')
  ax5.tick_params(which='both',direction='out')
  
  fig.tight_layout()
  mpl.rc('font',size=5)
  plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
  return




def sp_shape(pulses,events,store,obs):

  fig = plt.figure()

  mpl.rc('font',size=5)
 
  for i in range(0,pulses.shape[0]):
  
    puls = pulses.iloc[i]
    event = events.loc[events.Pulse==puls.name]

    sig = (event.Sigma/event.Sigma.max()*5)**4
  
    ax = plt.subplot2grid((2,5),(i/5,i%5))
    ax.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])  
    ax.errorbar(puls.Time_c, puls.DM_c, xerr=puls.dTime, yerr=puls.dDM, fmt='none', ecolor='r')
    
    ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(event.Sigma.max(),i))
    
  
  # Set common labels
  fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  fig.text(0.5, 0.95, obs, ha='center', va='center', fontsize=12)
  
  plt.savefig('{}'.format(store),format='png',bbox_inches='tight',dpi=200)
  
 
  return





def obs_top_candidates(top_candidates,color=True,size=True,store=False,incoherent=False): #top_candidates di tutti i beams
  
  if color:
    if incoherent:
      col_top = top_candidates.SAP
    else:
      col_top = top_candidates.BEAM
  else: 
    col_top = u'r' 
  
  if size: 
    sig_top = (top_candidates.Sigma/6.)**4
  else: sig=100.
    
  fig = plt.figure()
  plt.title("Best Candidates")
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
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
      bar.set_ticklabels(['{0:.0f}'.format(int(t)) for t in ticks])
      bar.ax.set_xlabel('beam',ha='left',labelpad=10)
      bar.update_ticks
      bar.ax.xaxis.set_ticks_position('top')
    
  ax1.plot([0,3600],[40.48,40.48],'k--')
  ax1.plot([0,3600],[141.68,141.68],'k--')
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,5,550])
  ax1.set_yscale('log')

  if not top_candidates.empty:
    if len(top_candidates.DM.unique())>1:
      ax2.hist(top_candidates.Sigma.tolist(),bins=100,histtype='step',color='k')
      ax2.set_xlabel('SNR')
      ax2.set_ylabel('Counts')
      ax2.set_yscale('log')
      
      ax3.hist(top_candidates.DM.tolist(),bins=300,histtype='stepfilled',color=u'k')
      ax3.set_xscale('log')
      ax3.set_xlabel('DM (pc/cm3)')
      ax3.set_ylabel('Counts')
      ax3.set_xlim(5,550)
    
    ax4.scatter(top_candidates.DM,top_candidates.Sigma,s=3,linewidths=[0.,],c=col_top)

    ax4.set_xscale('log')
    ax4.set_ylabel('SNR')
    ax4.set_xlabel('DM (pc/cm3)')
    ax4.axis([5,550,top_candidates.Sigma.min(),top_candidates.Sigma.max()+3.])
    ax4.plot([40.5,40.5],[0,top_candidates.Sigma.max()+3.],'k--')
    ax4.plot([141.68,141.68],[0,top_candidates.Sigma.max()+3.],'k--')
  
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
    
      
