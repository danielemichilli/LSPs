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
#import presto

import Utilities
from Parameters import *
import Paths
import RFIexcision



def alerts(pulses,folder):
  file_name = os.path.dirname(folder) + '/Candidates.log'
  cumSNR = pulses.groupby('Candidate').agg({'SAP': np.mean, 'BEAM': np.mean, 'DM': np.mean, 'Sigma': np.sum, 'N_events': np.size})
  cumSNR = cumSNR[cumSNR.Sigma>=10.]
  if not cumSNR.empty:
    cumSNR.to_csv(file_name,sep='\t',float_format='%.2f',columns=['SAP','BEAM','DM','N_events','Sigma'],header=['SAP','BEAM','DM','Num','Sigma'],index_label='Cand',mode='a')  
  return

def plot(((pulses,rfi,meta_data,events),folder,(sap,beam))):
  #pulses = gb[0]
  #rfi = gb[1]
  #meta_data = gb[2]
  #events = gb[3]
  #folder = gb[4]
  obs = os.path.basename(os.path.dirname(os.path.dirname(folder)))
  #sap = gb[5][0]
  #beam = gb[5][1]
  
  if pulses.empty: return

  alerts(pulses,folder)
  
  plt.clf()

  if beam == 12:
    sp_shape(pulses.head(10),events,'{}/SAP{}_BEAM{}/top_candidates(0-9).png'.format(folder,sap,beam),obs)
    sp_shape(pulses.iloc[10:20],events,'{}/SAP{}_BEAM{}/top_candidates(10-19).png'.format(folder,sap,beam),obs)
    sp_shape(pulses.iloc[20:30],events,'{}/SAP{}_BEAM{}/top_candidates(20-29).png'.format(folder,sap,beam),obs)
    plt.clf()
    sp_plot(pulses.iloc[30:],rfi,meta_data,pulses.head(30),sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
    
  else:
    sp_shape(pulses.head(10),events,'{}/SAP{}_BEAM{}/top_candidates.png'.format(folder,sap,beam),obs)
    plt.clf()
    sp_plot(pulses.iloc[10:],rfi,meta_data,pulses.head(10),sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
  
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

  #if not events.empty: 
    #ax1.scatter(events.Time, events.DM, facecolors='none', s=sig, c='k', linewidths=[0.5,])
      
  ax1.plot([0,3600],[40.48,40.48],'k--')
  ax1.plot([0,3600],[141.68,141.68],'k--')
    
  if not top_candidates.empty: ax1.scatter(top_candidates.Time, top_candidates.DM, s=sig_top, linewidths=[0.,], c=fill, marker='*')

  ax1.set_yscale('log')
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('DM (pc/cm3)')
  ax1.axis([0,3600,DM_MIN,550])
  
  mpl.rc('font', size=5)
  for i in range(0,top_candidates.shape[0]):
    if top_candidates.Time.iloc[i]>100:
      ax1.annotate(i,xy=(top_candidates.Time.iloc[i]-60,top_candidates.DM.iloc[i]),horizontalalignment='right',verticalalignment='center')
  hist = ax2.hist(pulses.DM.tolist()+top_candidates.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',range=(DM_MIN,550))
  ax2.set_xscale('log')
  ax2.set_xlabel('DM (pc/cm3)')
  ax2.set_ylabel('Counts')
  ax2.set_xlim(DM_MIN,550)
  ax2.plot([40.5,40.5],[0,hist[0].max()+10],'k--')
  ax2.plot([141.7,141.7],[0,hist[0].max()+10],'k--')
  
  ax3.scatter(pulses.DM,pulses.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=5,vmax=10)
  ax3.scatter(top_candidates.DM,top_candidates.Sigma,s=15.,linewidths=[0.,],c=fill,marker='*')
  ax3.set_xscale('log')
  ax3.set_ylabel('SNR')
  ax3.set_xlabel('DM (pc/cm3)')
  limit = max(pulses.Sigma.max(),top_candidates.Sigma.max())
  ax3.axis([DM_MIN,550,pulses.Sigma.min(),limit+3.])
  ax3.plot([40.5,40.5],[0,limit+3.],'k--')
  ax3.plot([141.7,141.7],[0,limit+3.],'k--')
  mpl.rc('font', size=3.5)
  for i in range(0,top_candidates.shape[0]):
    ax3.annotate(i,xy=(top_candidates.DM.iloc[i]/1.15,top_candidates.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center')
    
  #SNR = pulses.groupby('DM',sort=False).Sigma.sum()
  #ax4.plot(SNR.index,SNR,'k',linewidth=3.)
  hist = ax4.hist(pulses.DM.tolist()+top_candidates.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist()+top_candidates.Sigma.tolist(),range=(DM_MIN,550))
  ax4.set_xscale('log')
  ax4.set_xlabel('DM (pc/cm3)')
  ax4.set_ylabel('Cumulative SNR')
  ax4.set_xlim(DM_MIN,550)
  ax4.plot([40.5,40.5],[0,hist[0].max()+10],'k--')
  ax4.plot([141.7,141.7],[0,hist[0].max()+10],'k--')  

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
 
  for i,(idx,puls) in enumerate(pulses.iterrows()):
  
    event = events[events.Pulse==idx]

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
  ax1.axis([0,3600,DM_MIN,550])
  ax1.set_yscale('log')

  if not top_candidates.empty:  #PROBLEMA!!
    if len(top_candidates.DM.unique())>1:
      hist = ax2.hist(top_candidates.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',range=(DM_MIN,550))
      ax2.set_xscale('log')
      ax2.set_xlabel('DM (pc/cm3)')
      ax2.set_ylabel('Counts')
      ax2.set_xlim(DM_MIN,550)
      ax2.plot([40.5,40.5],[0,hist[0].max()+10],'k--')
      ax2.plot([141.7,141.7],[0,hist[0].max()+10],'k--')  
      
      ax3.scatter(top_candidates.DM,top_candidates.Sigma,s=3,linewidths=[0.,],c=col_top)
      ax3.set_xscale('log')
      ax3.set_ylabel('SNR')
      ax3.set_xlabel('DM (pc/cm3)')
      ax3.axis([DM_MIN,550,top_candidates.Sigma.min(),top_candidates.Sigma.max()+3.])
      ax3.plot([40.5,40.5],[0,top_candidates.Sigma.max()+3.],'k--')
      ax3.plot([141.7,141.7],[0,top_candidates.Sigma.max()+3.],'k--')

      hist = ax4.hist(top_candidates.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',weights=top_candidates.Sigma.tolist(),range=(DM_MIN,550))
      #ax4.plot(SNR.index,SNR,'k',linewidth=3.)
      ax4.set_xscale('log')
      ax4.set_xlabel('DM (pc/cm3)')
      ax4.set_ylabel('Cumulative SNR')
      ax4.set_xlim(DM_MIN,550)
      ax4.plot([40.5,40.5],[0,hist[0].max()+10],'k--')
      ax4.plot([141.7,141.7],[0,hist[0].max()+10],'k--')  
  
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
    





def DynamicSpectrum(pulses,idL,sap,beam,store):    #Creare una copia di pulses quando si chiama la funzione!
  if beam==12: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  filename = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(filename): return
  
  
  pulses.Sample[pulses.DM>141.71] *= 4
  pulses.Sample[(pulses.DM>40.47)&(pulses.DM<=141.71)] *= 2
  
  #controllare che questo vada dopo downsampling correction!
  header = Utilities.read_header(filename)
  MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
  v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
  bin_idx = np.int(np.round(1./v))
  pulses.Sample += pulses.Sample/bin_idx
  
  
  
  #fig = plt.figure()
  #mpl.rc('font',size=5)
  
  freq = np.arange(151,118,-1,dtype=np.float)
  offset = 300
  
  #for i,(idx,puls) in enumerate(pulses.iterrows()):
    #spectrum, bin_start, bin_end = Utilities.read_fits(filename,puls.DM,puls.Sample,offset)
    ##if RFIexcision.Multimoment() > FILTERS['Multimoment']: 
      ##pulses.Pulse += 1
      ##continue
    #extent = [bin_start*RES,bin_end*RES,119,151]
    #ax = plt.subplot2grid((2,5),(i/5,i%5))
    #ax.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
    #time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2)
    #ax.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
    ##ax.plot((time[0],time[0]+puls.Sample*RES),(freq[0],freq[0]),'r',(time[-1],time[-1]+puls.Sample*RES),(freq[-1],freq[-1]),'r')
    #ax.axis(extent)
    #ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(puls.Sigma,i))



  puls = pulses
  spectrum, bin_start, bin_end = Utilities.read_fits(filename,puls.DM,puls.Sample,offset)
  extent = [bin_start*RES,bin_end*RES,119,151]
  plt.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
  time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2) + bin_start*RES
  plt.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
  plt.axis(extent)
  

  
  # Set common labels
  #fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  #fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  #fig.text(0.5, 0.95, str(idL), ha='center', va='center', fontsize=12)
  
  #plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  
  return 


