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
import matplotlib.gridspec as gridspec
import logging
from matplotlib.backends.backend_pdf import PdfPages
import presto
import subprocess
import shutil

import Utilities
from Parameters import *
import Paths
import RFIexcision
from Paths import *

mpl.rc('font',size=5)



def output(idL,pulses,meta_data,candidates):
  pulses.sort('Sigma',ascending=False,inplace=True)
  candidates.sort(['Rank','Sigma'],ascending=False,inplace=True)
  
  plt.clf()
  out_dir = TMP_FOLDER.format(idL) + '/timeseries'
  if os.path.isdir(out_dir): shutil.rmtree(out_dir)
  for idx_c, cand in candidates.iterrows():
    sap = int(cand.SAP)
    beam = int(cand.BEAM)
    events = pd.read_hdf('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(idL)),'events',where=['(SAP == sap) & (BEAM == beam)'])
    store = '{}/sp/candidates/{}.pdf'.format(WRK_FOLDER.format(idL), cand.id)
    with PdfPages(store) as pdf:
      pulses_cand = pulses[pulses.Candidate == idx_c]
      beam_plot(pdf, cand, pulses_cand, pulses, meta_data)
      
      for i, (idx_p, puls) in enumerate(pulses_cand.head(5).iterrows()):
        ev = events[events.Pulse == idx_p]
        puls_plot(pdf, puls, ev, idL, i)
    
    plt.close('all')
    if os.path.isdir(out_dir): shutil.rmtree(out_dir)
    
  return



def beam_plot(pdf, cand, pulses, pulses_all, meta_data):
  sap = int(cand.SAP)
  beam = int(cand.BEAM)
  pulses_beam = pulses_all[(pulses_all.SAP==sap) & (pulses_all.BEAM==beam)]
    
  ax1 = plt.subplot2grid((3,6),(0,0), rowspan=2)
  ax2 = plt.subplot2grid((3,6),(0,1), colspan=5, rowspan=2)
  ax3 = plt.subplot2grid((3,6),(2,0), colspan=2)
  ax4 = plt.subplot2grid((3,6),(2,2), colspan=2)
  ax5 = plt.subplot2grid((3,6),(2,4), colspan=2)

  meta_data_plot(ax1,meta_data[(meta_data.SAP==sap)&(meta_data.BEAM==beam)],pulses,cand)
  scatter_beam(ax2, pulses, pulses_beam, cand)

  if not pulses_beam.empty: scatter_SNR(ax4,pulses,pulses_beam,cand)
  try: 
    hist_DM(ax3,pulses_beam,cand)
    hist_SNR(ax5,pulses_beam,cand)
  except ValueError: pass
    
  plt.tight_layout()
  pdf.savefig(bbox_inches='tight',dpi=200)
  return



def puls_plot(pdf, puls, ev, idL, i):
  ax1 = plt.subplot2grid((2,6),(0,0))
  ax5 = plt.subplot2grid((2,6),(0,1))
  ax6 = plt.subplot2grid((2,6),(0,2), colspan=2)
  ax7 = plt.subplot2grid((2,6),(0,4), colspan=2)
  ax2 = plt.subplot2grid((2,6),(1,0))
  ax3 = plt.subplot2grid((2,6),(1,1), colspan=2)
  ax4 = plt.subplot2grid((2,6),(1,3), colspan=3)
  
  puls_meta_data(ax1, puls, ev.Pulse.iloc[0], i)
  puls_DM_Time(ax2, ev, puls)
  puls_SNR_DM(ax3, ev)
  if puls.BEAM > 12: puls_heatmap(ax4, puls, idL)
  else: plot_not_valid(ax4)
  flag = puls_dynSpec(ax5, ax6, puls, idL)
  if flag == -1:
    plot_not_valid(ax5)
    plot_not_valid(ax6)
  dir_ts = TMP_FOLDER.format(idL) + '/timeseries/{}_SAP{:.0f}_BEAM{:.0f}_DM{{0:.2f}}.dat'.format(idL, puls.SAP, puls.BEAM)
  flag = puls_dedispersed(ax7, puls, dir_ts, idL)
  if flag == -1:
    plot_not_valid(ax7)

  plt.tight_layout()
  pdf.savefig(bbox_inches='tight',dpi=200)
  return



def plot_not_valid(ax):
  #ax.text(.5,.5,'Plot not valid', size=50., horizontalalignment='center', verticalalignment='center', clip_on=True)
  ax.plot([0,1], [0,1], 'k-', lw=2.)
  ax.plot([0,1], [1,0], 'k-', lw=2.)
  ax.set_xticks([])
  ax.set_yticks([])
  return



def scatter_beam(ax, pulses, pulses_beam, cand):
  def circle_size(values):
    new_val = np.clip(values,6.5,20)
    m = 31.85
    q = -137.025
    return new_val * m + q
  sig = circle_size(pulses_beam.Sigma)

  colors = ['g','y','r']
  col = pulses_beam.Pulse.replace([0,1,2], colors)

  ax.scatter(pulses_beam.Time, pulses_beam.DM, c=col, s=sig, edgecolor='w', lw=.2, zorder=2)
  
  def scatter_legend(colors):
    scatter_list = []
    for c in colors:
      scatter_list.append( ax.scatter([],[],linewidths=[0,],c=c))
    return scatter_list
  ax.legend(scatter_legend(colors), ('Best','Medium','Worst'), loc='upper right')
  ax.axhline(cand.DM, c='dodgerblue', ls='-', zorder=1)
  
  ax.axhline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axhline(141.68,c='k',ls='--',lw=.1, zorder=1)
  ax.set_yscale('log')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cc)')
  ax.axis([0,3600,DM_MIN,550])
  
  sig = circle_size(pulses.Sigma)
  ax.scatter(pulses.Time, pulses.DM, s=sig, linewidths=[0.,], marker='*', c='w', zorder=3)
  
  top = pulses.iloc[:10]
  for i in range(top.shape[0]):
    ax.annotate(str(i),xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',verticalalignment='center',color='dodgerblue',size=5, zorder=4)
  ax.tick_params(which='both',direction='out')
  
  return



def meta_data_plot(ax,meta_data,pulses,cand):
  ax.axis([0,10,0,10])
  ax.annotate('File: {}'.format(meta_data.File.iloc[0]), xy=(0,9))
  ax.annotate('RA, DEC: {0:.8}, {1:.8}'.format(meta_data.RA.iloc[0],meta_data.DEC.iloc[0]), xy=(0,8))
  ax.annotate('Epoch (MJD): {0:.11}'.format(meta_data.Epoch.iloc[0]), xy=(0,7))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(cand.DM), xy=(0,6))
  ax.annotate('Sigma (cum.): {0:.1f}'.format(cand.Sigma), xy=(0,5))
  if cand.N_pulses == 1:
    ax.annotate('Time (s): {0:.2f}'.format(pulses.Time.iloc[0]), xy=(0,3))
    ax.annotate('Duration (ms): {0:.0f}'.format(pulses.Duration.iloc[0]*1000), xy=(0,2))
    ax.annotate('N. events: {}'.format(pulses.N_events.iloc[0]), xy=(0,1))
  else:
    ax.annotate('dDM (pc/cm2): {0:.2f}'.format(pulses.DM.max()-pulses.DM.min()), xy=(0,4))
    ax.annotate('Period (s): {0:.3f}'.format(cand.Period), xy=(0,3))
    ax.annotate('Period err. (s): {0:.3f}'.format(cand.Period_err), xy=(0,2))
    ax.annotate('N. pulses: {0:.0f}'.format(cand.N_pulses), xy=(0,1))
  ax.axis('off')
  return



def hist_DM(ax,pulses,cand):
  ax.axvline(cand.DM, c='dodgerblue', ls='-', linewidth=.2, zorder=3)
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',range=(DM_MIN,550), zorder=2)
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Counts')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(141.68,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def hist_SNR(ax,pulses,cand):
  ax.axvline(cand.DM, c='dodgerblue', ls='-', linewidth=.2, zorder=3)
  ax.hist(pulses.DM.tolist(),bins=550-DM_MIN,histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist(),range=(DM_MIN,550), zorder=2)
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Cumulative SNR')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(141.68,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def scatter_SNR(ax, pulses, pulses_beam,cand):
  ax.axvline(cand.DM, c='dodgerblue', ls='-', linewidth=.2, zorder=2)
  colors = ['g','y','r']
  col = pulses_beam.Pulse.replace([0,1,2], colors)  
  ax.scatter(pulses_beam.DM,pulses_beam.Sigma,c=col,s=20.,linewidths=[0.,], zorder=3)
  ax.set_xscale('log')
  ax.set_ylabel('SNR')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_xlim((3., 550.))
  ax.set_ylim((6.5, pulses_beam.Sigma.max()+1))
  ax.axvline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(141.68,c='k',ls='--',lw=.1, zorder=1)
  top = pulses.iloc[:10]
  for i in range(0,top.shape[0]):
    ax.annotate(i,xy=(top.DM.iloc[i]/1.15,top.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center',size=4,weight='medium')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def puls_DM_Time(ax, event, puls):
  #def circle_size(val):
    #m = 4824.
    #q = 17.59
    #return val * m + q
  #sig = circle_size(event.Duration)
  #ax.scatter(event.Time, event.DM, facecolors='none', s=sig, c='k', linewidths=[0.5,], zorder=2)  
  ax.scatter(event.Time, event.DM, s=20., marker='o', c='k', linewidths=[0.,], zorder=1)  
  #ax.errorbar(puls.Time, puls.DM, xerr=puls.Duration/2, yerr=puls.dDM/2, lw=.2, fmt='none', ecolor='r', zorder=2)
  ax.scatter(puls.Time, puls.DM, lw=.2, marker='x', color='r', zorder=2)
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cm3)')
  return



def puls_SNR_DM(ax, event):
  ax.plot(event.DM, event.Sigma, 'k')
  ax.set_xlabel('DM (pc/cm3)')    
  ax.set_ylabel('SNR')
  
  ax2 = ax.twinx()
  ax2.plot(event.DM, event.Duration*1000, 'r')
  ax2.set_ylabel('Duration (ms)', color='r')
  for tl in ax2.get_yticklabels():
    tl.set_color('r')
  return



def puls_heatmap(ax, puls, idL, pulseN=False):
  ra = np.array([ 
        499,  499,  622,  621,  499,  377,  376,  499,  624,  748,  744,
        742,  620,  499,  379,  257,  254,  251,  375,  499,  625,  750,
        873,  869,  865,  862,  740,  619,  499,  380,  259,  137,  134,
        129,  126,  249,  374,  499,  627,  752,  876, 1000,  995,  990,
        985,  981,  859,  738,  618,  499,  381,  261,  140,   18,   14,
          9,    4,    0,  123,  247,  372])
  dec = np.array([ 
        500,  625,  562,  437,  375,  437,  562,  750,  687,  625,  500,
        375,  312,  250,  312,  375,  500,  625,  687,  875,  812,  750,
        687,  562,  437,  312,  250,  187,  125,  187,  250,  312,  437,
        562,  687,  750,  812, 1000,  937,  875,  812,  750,  625,  500,
        375,  250,  187,  125,   62,    0,   62,  125,  187,  250,  375,
        500,  625,  750,  812,  875,  937])
  
  dDM = 0.4
  dm_l = float(puls.DM - dDM/2. - 0.001)
  dm_h = float(puls.DM + dDM/2. + 0.001)
  t_l = float(puls.Time - 2. * puls.Duration)
  t_h = float(puls.Time + 2. * puls.Duration)
  sap = int(puls.SAP)
  select = '(SAP == sap) and ((DM > dm_l) and (DM < dm_h)) and ((Time >= t_l) and (Time <= t_h))'
  try:
    events = pd.read_hdf('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(idL)), 'events', where=select)
  except ValueError:
    events = pd.read_hdf('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(idL)), 'events')
    events = events[(events.SAP == sap) & ((events.DM > dm_l) & (events.DM < dm_h)) & ((events.Time >= t_l) & (events.Time <= t_h))]
  SNR = events.groupby('BEAM').Sigma.sum()
  ind = pd.Series(np.zeros(61))
  ind.index += 13
  SNR = SNR.reindex_like(ind)

  plot = ax.scatter(ra,dec,s=400,edgecolor='none',c=SNR,cmap='hot_r')
  bar = plt.colorbar(plot, ax=ax)
  
  bar.set_label('Cumulative SNR')
  ax.set_xlabel('RA (rel.)')
  ax.set_ylabel('DEC (rel.)')
  if pulseN: 
    ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN))
    ax.annotate('DM: {:.2f}$\pm$0.2, Time: {:.2f}$\pm${:.2f}'.format(puls.DM,puls.Time,puls.Duration), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  else:
    ax.annotate('DM: {:.2f} - {:.2f}, Time: {:.2f} - {:.2f}'.format(dm_l,dm_h,t_l,t_h), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  
  beam = puls.BEAM
  ax.scatter(ra[beam-13],dec[beam-13],s=300,linewidths=[0.,],marker='*',c='w')
  [ax.annotate(str(i+13),(ra[i],dec[i]), horizontalalignment='center', verticalalignment='center', color='m') for i in range(0,61)]

  ax.set_xlim(-100,1100)
  ax.set_ylim(-100,1100)
  return



def puls_meta_data(ax, puls, idx, i):
  ax.axis([0,10,0,10])
  ax.annotate("Pulse n. {}".format(i), xy=(0,9)) 
  ax.annotate('Pulse code: {}'.format(idx), xy=(0,8))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(puls.DM), xy=(0,7))
  ax.annotate('dDM (pc/cm2): {0:.2f}'.format(puls.dDM), xy=(0,6))
  ax.annotate('Time (s): {0:.2f}'.format(puls.Time), xy=(0,5))
  ax.annotate('Sigma: {0:.1f}'.format(puls.Sigma), xy=(0,4))
  ax.annotate('Duration (ms): {0:.0f}'.format(puls.Duration*1000), xy=(0,3))
  ax.annotate('N events: {0:.0f}'.format(puls.N_events), xy=(0,2))
  ax.annotate('Rank: {0:.0f}'.format(puls.Pulse), xy=(0,1))

  ax.axis('off')
  return



def puls_dynSpec(ax1, ax2, puls, idL):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam==12: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  filename = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(filename): return -1

  if puls.DM>141.71: sample = puls.Sample * 4
  elif puls.DM>40.47: sample = puls.Sample * 2
  else: sample = puls.Sample

  filetype, header = Utilities.read_header(filename)
  MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
  v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
  sample += np.round(sample*v).astype(int)

  duration = np.int(np.round(puls.Duration/RES))
  spectra_border = 20
  offset = duration*spectra_border

  #Load the spectrum
  mask_name = MASK_FILE.format(idL=idL,sap=sap,beam=beam)
  spectrum = Utilities.read_fits(filename,puls.DM.copy(),sample.copy(),duration,offset,RFI_reduct=True,mask=mask_name)
  extent = [(sample-offset)*RES,(sample+duration+offset)*RES,F_MIN,F_MAX]

  if len(spectrum) == 0: return -1
  
  #Probabilmente extent da rifare per spettro dispersed
  ax2.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
  #Plottare contorni!
  ax2.axis(extent)
  ax2.set_xlabel('Time (s)')
  ax2.set_ylabel('Frequency (MHz)')  

  def dedispersion(spectrum):
    freq = np.linspace(F_MIN,F_MAX,2592)
    time = (4149 * puls.DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
    for i in range(spectrum.shape[1]):
      spectrum[:,i] = np.roll(spectrum[:,i], time[i])
    spectrum = spectrum[:2*offset+duration]
  
    spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0]/duration,duration,spectrum.shape[1]]),axis=1)
    spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0],spectrum.shape[1]/32,32]),axis=2)
    return spectrum
  
  spectrum = dedispersion(spectrum)
     
  ax1.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
  ax1.scatter((sample+duration/2)*RES,F_MIN+1,marker='^',s=100,c='r',lw=0.)
  ax1.axis(extent)
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('Frequency (MHz)')
      
  return 0



def load_ts(puls, filename, idL):
  k = 4148.808 * (F_MAX**-2 - F_MIN**-2) / RES  #Theoretical factor between time and DM

  bin_peak = int(puls['Sample'])
  DM_peak = puls['DM']
  if DM_peak < 40.52: 
    duration = int(np.round(puls['Duration'] / RES))
    DM_res = 0.01
  elif DM_peak < 141.77: 
    duration = int(np.round(puls['Duration'] / RES / 2.))
    DM_res = 0.05
    k /= 2.
  else: 
    duration = int(np.round(puls['Duration'] / RES / 4.))
    DM_res = 0.1
    k /= 4.

  nDMs = int(np.round(puls['dDM'] / DM_res * 2.))  #Valutare se estendere oltre larghezza puls
  DM_range = DM_peak - (nDMs/2 - np.arange(nDMs)) * DM_res 
  
  nProfBins = 3
  scrunch_fact = int(np.round(duration / float(nProfBins)))
  if scrunch_fact < 1: scrunch_fact = 1  
  
  nPlotBins = int(np.ceil((DM_range[0] - DM_range[-1]) * k / duration * 1.5 ))
  nBins = nPlotBins * nProfBins
  data = np.zeros((nDMs,nBins))

  bin_start = bin_peak - nBins/2 * scrunch_fact

  out_dir = TMP_FOLDER.format(idL) + '/timeseries'
  FNULL = open(os.devnull, 'w')
  for j,DM in enumerate(DM_range):
    if not os.path.isfile(filename.format(DM)):
      error = subprocess.call(['sh', '/home/danielem/spdspsr_DANIELE.sh', idL, str(sap), str(beam), '{:.2f}'.format(DM), '1', out_dir], stdout=FNULL, stderr=FNULL)
    
    try:
      ts = np.memmap(filename.format(DM), dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
    except IOError: ts = np.zeros(nBins) + np.nan

    data[j] = ts
    
  params = {'bins_out': nBins*scrunch_fact, 'bin_start': bin_start, 'scrunch_fact': scrunch_fact, 'duration': duration, 'DM_range': DM_range, 'k': k}
  
  return data, params



def puls_dedispersed(ax, puls, filename, idL, pulseN=False):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam==12: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  raw_dir = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(raw_dir): return -1
    
  data, params = load_ts(puls, filename, idL)
  
  #Image plot
  ax.imshow(data,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=[0,params['bins_out'],params['DM_range'][0], params['DM_range'][-1]])
  ax.set_ylim((params['DM_range'][0], params['DM_range'][-1]))
  ax.set_ylabel('DM (pc/cc)')
  if pulseN: ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN), y=1.08)

  #Time axis
  ax2 = ax.twiny()
  x = np.arange(params['bins_out'])
  ax2.set_xlim((0,params['bins_out']))
  xticks_pos = x[::(params['bins_out'])/6]
  xticks_val = (xticks_pos + params['bin_start']) * RES
  ax2.set_xticks(xticks_pos, ["%.3f" % n for n in xticks_val])
  ax2.set_xlabel('Time (s)')
  
  #Plot contours
  y = (x - params['bins_out'] / 2.) / params['k'] + puls.DM
  ax2.plot(x, y,'r--')
  ax2.axvline((params['bins_out'] - params['scrunch_fact']) / 2., color='r', ls='--')

  #Sample axis
  ax2.set_xlim((0,params['bins_out']))
  ax2.set_xticks(xticks_pos, ((xticks_pos + params['bin_start']) % 10000).tolist() )
  ax2.set_xlabel('Sample - {}0000'.format(params['bin_start'] / 10000))
  
  #Inset profile
  def inset(filename, params, puls):
    nPlotBins = 20
    nProfBins = 5
    nBins = nPlotBins * nProfBins
    scrunch_fact = int(np.round(params['duration'] / float(nProfBins)))
    if scrunch_fact < 1: scrunch_fact = 1
    x = np.arange(nBins)
    bin_start = int(puls.Sample) - nBins/2 * scrunch_fact
  
    try:
      ts = np.memmap(filename.format(puls.DM), dtype=np.float32, mode='r', offset=bin_start*4, shape=(params['bins_out'],))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
    except IOError: x = ts = None
    return x, ts
  
  x, ts = inset(filename, params, puls)
  ax.axis([.7, .7, .15, .15])
  ax.plot(x, ts, 'k')
  ax.set_xticks([])
  ax.set_yticks([])

  return 0
