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

mpl.rc('font',size=5)


def alerts(pulses,folder):
  file_name = os.path.dirname(folder) + '/Candidates.log'
  cumSNR = pulses.groupby('Candidate').agg({'SAP': np.mean, 'BEAM': np.mean, 'DM': np.mean, 'Sigma': np.sum, 'N_events': np.size})
  cumSNR = cumSNR[cumSNR.Sigma>=10.]
  if not cumSNR.empty:
    cumSNR.to_csv(file_name,sep='\t',float_format='%.2f',columns=['SAP','BEAM','DM','N_events','Sigma'],header=['SAP','BEAM','DM','Num','Sigma'],index_label='Cand',mode='a')  
  return



def sp_plot(pulses,rfi,meta_data,sap,beam,store):
  plt.clf()

  cmap = plt.get_cmap('gist_heat_r')
  fill = u'b'
  square = u'g'
    
  fig = plt.figure()
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  ax5 = plt.subplot2grid((3,4),(0,3))

  scatter_beam(ax1,pulses,cmap,rfi=rfi)
  try: hist_DM(ax2,pulses)
  except ValueError: pass
  scatter_SNR(ax3,pulses,cmap)
  try: hist_SNR(ax4,pulses)
  except ValueError: pass
  meta_data_plot(ax5,meta_data)
  
  fig.tight_layout()
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  return


def sp_shape(pulses,events,store,obs):
  plt.clf()
  fig = plt.figure()
  
  for i,(idx,puls) in enumerate(pulses.iterrows()):
    event = events[events.Pulse==idx]
    ax = plt.subplot2grid((2,5),(i/5,i%5))
    puls_DM_Time(ax,event,puls)
    ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(puls.Sigma,i))
    
  # Set common labels
  fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  fig.text(0.5, 0.95, '{} - SAP{}_BEAM{}'.format(obs,pulses.SAP.unique()[0],pulses.BEAM.unique()[0]), ha='center', va='center', fontsize=12)
  
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  return


def obs_top_candidates(top_candidates,store,incoherent=False):
  plt.clf()
  
  if incoherent:
    col = top_candidates.SAP
    num = top_candidates.SAP.unique().size
  else:
    col = top_candidates.BEAM
    num = top_candidates.BEAM.unique().size
  
  if num > 1: cmap = discrete_cmap(num, 'spectral')
  else: cmap = plt.get_cmap('gist_heat_r')
  fig = plt.figure()
  plt.title("Best Candidates")
  
  ax1 = plt.subplot2grid((3,4),(1,0),colspan=4,rowspan=2)
  ax2 = plt.subplot2grid((3,4),(0,0))
  ax3 = plt.subplot2grid((3,4),(0,1))
  ax4 = plt.subplot2grid((3,4),(0,2))
  
  scatter_beam(ax1,top_candidates,cmap,col=col,legend=num)
  
  try: hist_DM(ax2,top_candidates)
  except ValueError: pass
  scatter_SNR(ax3,top_candidates,cmap,col=col,with_legend=True)
  try: hist_SNR(ax4,top_candidates)
  except ValueError: pass
  
  fig.tight_layout()
  plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  
  return





#DA RIMUOVERE (o spostare in utilities)
def DynamicSpectrum(pulses,idL,sap,beam,store):    #Creare una copia di pulses quando si chiama la funzione!
  plt.clf()
  
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
  
  vmin,vmax = Utilities.color_range(spectrum)
  
  plt.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent,vmin=vmin,vmax=vmax)     #provare anche pcolormesh
  time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2) + bin_start*RES
  plt.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
  plt.axis(extent)
  
  
  
  
  
  #allineamento perfetto
  #problema di contrasto: provare a diminuire range di colori, azzerare minori di soglia, mediare in duration del pulse, etc
  #prova di de-dispersione (forse inutile):
  freq = np.linspace(F_MIN,F_MAX,2430)
  time = (4149 * DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
  for i in range(a.shape[1]):
    a[:,i] = np.roll(a[:,i], time[i])

  
  
  
  
  # Set common labels
  #fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  #fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  #fig.text(0.5, 0.95, str(idL), ha='center', va='center', fontsize=12)
  
  #plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  
  return 





def single_candidates(events,pulses,cands,meta_data,idL,folder):
  for i,(idx,cand) in enumerate(cands.iterrows()):
    #print cand
    #print pulses
    puls = pulses[pulses.Candidate==idx]
    #print puls
    #print events
    event = events[events.Pulse==puls.index[0]]
    sap = puls.SAP.iloc[0]
    beam = puls.BEAM.iloc[0]
    
    plt.clf()
    fig = plt.figure()
    
    ax1 = plt.subplot2grid((2,3),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,3),(0,2))
    ax3 = plt.subplot2grid((2,3),(1,0))
    ax4 = plt.subplot2grid((2,3),(1,1))
    ax5 = plt.subplot2grid((2,3),(1,2))
    
    scatter_beam(ax1,pulses,'gist_heat_r')
    ax1.scatter(puls.Time, puls.DM, s=100, linewidths=[0.,], marker='*')
    meta_data_puls(ax2,meta_data[(meta_data.SAP==sap)&(meta_data.BEAM==beam)],puls)
    ax3.plot(event.DM, event.SNR, 'k')
    puls_DM_Time(ax4,event,puls)
    DynamicSpectrum(ax5,pulses.copy(),idL,sap,beam)
    
    fig.tight_layout()
    plt.savefig('{}{}/sp/files/prova_{}.png'.format(folder,idL,i),format='png',bbox_inches='tight',dpi=200)

  return



#FINIRE!!
def repeated_candidates(events,pulses,cands,idL,folder):
  for i,(idx,cand) in enumerate(cands):
    puls = pulses[pulses.Candidate==cand.index]
    event = events[events.Pulse==puls.index[0]]
    
    plt.clf()
    fig = plt.figure()
    
    ax1 = plt.subplot2grid((2,3),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,3),(0,2))
    ax3 = plt.subplot2grid((2,3),(1,0))
    ax4 = plt.subplot2grid((2,3),(1,1))
    ax5 = plt.subplot2grid((2,3),(1,2))
    
    scatter_beam(ax1,pulses,'gist_heat_r')
    ax1.axhline(puls.DM,ls='--')
    meta_data_puls(ax2,meta_data,puls)
    ax3.plot(event.DM, event.SNR, 'k')
    puls_DM_Time(ax4,event,puls)
    DynamicSpectrum(ax5,pulses.copy(),idL,puls.SAP.iloc[0],puls.BEAM.iloc[0])
    
    fig.tight_layout()
    plt.savefig('{}{}/sp/files/prova_{}.png'.format(folder,idL,i),format='png',bbox_inches='tight',dpi=200)

  return






def scatter_beam(ax,pulses,cmap,col=None,rfi=False,legend=False):
  sig = np.clip(np.log(pulses.Sigma-5.5)*400+100,100,1200)
  if isinstance(col,type(None)): col = sig

  if legend:
    main_plt = ax.scatter(pulses.Time, pulses.DM, c=col, s=sig, cmap=cmap, linewidths=[0.,])
    ticks = np.linspace(col.min(),col.max(),num=legend)
    bar = plt.colorbar(main_plt,ticks=ticks,ax=ax)
    bar.set_ticklabels(['{0:.0f}'.format(int(t)) for t in ticks])
    bar.ax.set_xlabel('sap',ha='left',labelpad=10)
    bar.update_ticks
    bar.ax.xaxis.set_ticks_position('top')
  else: ax.scatter(pulses.Time, pulses.DM, c=col, s=sig, cmap=cmap, linewidths=[0.,], vmin=SNR_MIN, vmax=15)
  if isinstance(rfi,pd.DataFrame): ax.scatter(rfi.Time, rfi.DM, s=5., c=u'k', marker='+', linewidths=[0.4,])
  
  ax.axhline(40.48,c='k',ls='--',lw=.1)
  ax.axhline(141.68,c='k',ls='--',lw=.1)
  ax.set_yscale('log')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cm3)')
  ax.axis([0,3600,DM_MIN,550])
  
  top = pulses.iloc[:10]
  for i in range(top.shape[0]):
    ax.annotate(i,xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',verticalalignment='center',color='dodgerblue',size=7,weight='bold')
  ax.tick_params(which='both',direction='out')
  
  return


def hist_DM(ax,pulses):
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',range=(DM_MIN,550))
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Counts')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def scatter_SNR(ax,pulses,cmap,col=None,with_legend=False):
  if isinstance(col,type(None)): col = np.clip(np.log(pulses.Sigma-5.5)*400+100,100,1200)
  if with_legend: ax.scatter(pulses.DM,pulses.Sigma,c=col,s=6.,cmap=cmap,linewidths=[0.,])
  else: ax.scatter(pulses.DM,pulses.Sigma,c=col,s=3.,cmap=cmap,linewidths=[0.,],vmin=SNR_MIN,vmax=15)
  ax.set_xscale('log')
  ax.set_ylabel('SNR')
  ax.set_xlabel('DM (pc/cm3)')
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  top = pulses.iloc[:10]
  for i in range(0,top.shape[0]):
    ax.annotate(i,xy=(top.DM.iloc[i]/1.15,top.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center',size=4,weight='medium')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def hist_SNR(ax,pulses):
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist(),range=(DM_MIN,550))
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Cumulative SNR')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1)
  ax.axvline(141.68,c='k',ls='--',lw=.1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return


def meta_data_plot(ax,meta_data):
  ax.axis([0,10,0,7])
  ax.annotate('File: '+meta_data.File.iloc[0], xy=(0,6))
  ax.annotate('Telescope: '+meta_data.Telescope.iloc[0], xy=(0,5))
  ax.annotate('Instrument: '+meta_data.Instrument.iloc[0], xy=(0,4))
  ax.annotate('RA: '+meta_data.RA.iloc[0], xy=(0,3))
  ax.annotate('DEC: '+meta_data.DEC.iloc[0], xy=(0,2))
  ax.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,1))
  ax.axis('off')
  return

def meta_data_puls(ax,meta_data,puls):
  ax.axis([0,10,0,9])
  ax.annotate('File: '+meta_data.File.iloc[0], xy=(0,8))
  ax.annotate('RA, DEC: '+meta_data.RA.iloc[0]+', '+meta_data.DEC.iloc[0], xy=(0,7))
  ax.annotate('Epoch (MJD): '+meta_data.Epoch.iloc[0], xy=(0,6))
  ax.annotate('DM (pc/cm2): '+puls.DM.iloc[0].astype(str), xy=(0,5))
  ax.annotate('Time (s): '+puls.Time.iloc[0].astype(str), xy=(0,4))
  ax.annotate('Sigma: '+puls.Sigma.iloc[0].astype(str), xy=(0,3))
  ax.annotate('Duration (ms): '+str(puls.Sigma.iloc[0]*1000), xy=(0,2))
  ax.annotate('N_events: '+puls.N_events.iloc[0].astype(str), xy=(0,1))
  ax.axis('off')
  return

def puls_DM_Time(ax,event,puls):
  #sig = (event.Sigma/event.Sigma.max()*5)**4
  #sig = np.clip(np.log(event.Sigma-5.5)*400+100,100,1200)
  #sig = event.Sigma/event.Sigma.max()*1000
  sig = np.clip(event.Sigma/event.Sigma.max()*1427-427,1,1000)
  ax.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k',linewidths=[0.5,])  
  ax.errorbar(puls.Time, puls.DM, xerr=puls.dTime/2, yerr=puls.dDM/2, fmt='none', ecolor='r')
  return


def discrete_cmap(N, base_cmap):
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)  
  
  

def DynamicSpectrum(ax,pulses,idL,sap,beam):
  plt.clf()
  
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
  
  freq = np.arange(151,118,-1,dtype=np.float)
  offset = 300

  puls = pulses
  spectrum, bin_start, bin_end = Utilities.read_fits(filename,puls.DM,puls.Sample,offset)
  extent = [bin_start*RES,bin_end*RES,119,151]
  
  vmin,vmax = Utilities.color_range(spectrum)
  
  ax.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent,vmin=vmin,vmax=vmax)
  time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2) + bin_start*RES
  ax.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
  ax.axis(extent)
  
  return 
