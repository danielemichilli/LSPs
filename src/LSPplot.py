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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import presto
import subprocess
import shutil
import waterfaller
import psrfits

import Utilities
from Parameters import *
import RFIexcision
from Paths import *

mpl.rc('font',size=5)



def output(idL, pulses, meta_data, candidates, db, inc=12):
  pulses.sort('Sigma',ascending=False,inplace=True)
  candidates.sort('Sigma',ascending=False,inplace=True)
  
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
      beam_plot(pdf, cand, pulses_cand, pulses, meta_data, events)
      
      for i, (idx_p, puls) in enumerate(pulses_cand.head(5).iterrows()):
        puls_plot(pdf, puls, events, idL, db, i, inc=inc)
    
    plt.close('all')
    if os.path.isdir(out_dir): shutil.rmtree(out_dir)
    
  return



def beam_plot(pdf, cand, pulses, pulses_all, meta_data, events):
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

  if not pulses_beam.empty: 
    scatter_SNR(ax4,pulses,events[events.Pulse.isin(pulses_beam.index)],cand)
  try: 
    hist_DM(ax3,pulses_beam,cand)
    hist_SNR(ax5,pulses_beam,cand)
  except ValueError: pass
    
  plt.tight_layout()
  pdf.savefig(bbox_inches='tight',dpi=200)
  return



def puls_plot(pdf, puls, events, idL, db, i, inc=12):
  ax1 = plt.subplot2grid((2,6),(0,0))
  ax5 = plt.subplot2grid((2,6),(0,1))
  ax6 = plt.subplot2grid((2,6),(0,2), colspan=2)
  ax7 = plt.subplot2grid((2,6),(0,4), colspan=2)
  ax2 = plt.subplot2grid((2,6),(1,0), colspan=2)
  ax3 = plt.subplot2grid((2,6),(1,2), colspan=2)
  ax4 = plt.subplot2grid((2,6),(1,4), colspan=2)

  ev = events[events.Pulse == puls.name]
  puls_meta_data(ax1, puls, ev.Pulse.iloc[0], i)
  puls_DM_Time(ax2, ev, events, puls)
  puls_SNR_DM(ax3, ev)
  puls_heatmap(ax4, puls, idL, db, inc=inc)
  flag = puls_dynSpec(ax5, ax6, puls, idL, inc=inc)
  if flag == -1:
    plot_not_valid(ax5)
    plot_not_valid(ax6)
  flag = puls_dedispersed(ax7, puls, idL, inc=inc)
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

  ax.scatter(pulses_beam.Time, pulses_beam.DM, c='k', s=sig, edgecolor='w', lw=.2, zorder=2)

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
    ax.annotate(str(i),xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',verticalalignment='center',color='dodgerblue',size=5,fontweight='bold',zorder=4)
  ax.tick_params(which='both',direction='out')
  
  return



def meta_data_plot(ax,meta_data,pulses,cand):
  ax.axis([0,10,0,10])
  ax.annotate('Cand ID: {}'.format(cand.name), xy=(0,10))
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
  ax.hist(pulses.DM.tolist(),bins=int(550-DM_MIN),histtype='stepfilled',color=u'k',range=(DM_MIN,550), zorder=2)
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
  ax.hist(pulses.DM.tolist(),bins=int(550-DM_MIN),histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist(),range=(DM_MIN,550), zorder=2)
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_ylabel('Cumulative SNR')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(141.68,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def scatter_SNR(ax, pulses, pulses_beam, cand):
  ax.axvline(cand.DM, c='dodgerblue', ls='-', linewidth=.2, zorder=2)
  ax.scatter(pulses_beam.DM,pulses_beam.Sigma,c='k',s=10.,linewidths=[0.,], zorder=3)
  ax.set_xscale('log')
  ax.set_ylabel('SNR')
  ax.set_xlabel('DM (pc/cm3)')
  ax.set_xlim((3., 550.))
  ax.set_ylim((6.5, pulses_beam.Sigma.max()+1))
  ax.axvline(40.48,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(141.68,c='k',ls='--',lw=.1, zorder=1)
  #top = pulses.iloc[:10]
  #for i in range(0,top.shape[0]):
  #  ax.annotate(i,xy=(top.DM.iloc[i]/1.15,top.Sigma.iloc[i]),horizontalalignment='right',verticalalignment='center',size=4,weight='medium')
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def puls_DM_Time(ax, event, all_events, puls):
  def circle_size(values):
    new_val = np.clip(values,6.5,20)
    m = 31.85
    q = -137.025
    return new_val * m + q
  sig = circle_size(event.Sigma)
  ax.scatter(event.Time, event.DM, s=sig, marker='o', c='k', linewidths=[0.1,], facecolors='none', zorder=1)
  #ax.scatter(event.Time, event.DM, s=20., marker='o', c='k', linewidths=[0.,], zorder=1)
  #ax.errorbar(puls.Time, puls.DM, xerr=puls.Duration/2, yerr=puls.dDM/2, lw=.2, fmt='none', ecolor='r', zorder=2)
  ax.scatter(puls.Time, puls.DM, marker='x', color='r', zorder=2)
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cm3)')
  
  DM_min = np.clip(puls.DM - 5., 3., 550.)
  DM_max = np.clip(puls.DM + 5., 3., 550.)
  k = 4148.808 #s-1
  delay = k * (F_MIN**-2 - F_MAX**-2)
  dDM = DM_max - DM_min
  dt = delay * dDM / 4.
  Time_min = puls.Time - dt
  Time_max = puls.Time + dt
  events = all_events[(all_events.DM > DM_min) & (all_events.DM < DM_max) & (all_events.Time > Time_min) & (all_events.Time < Time_max)]
  ax.scatter(events.Time, events.DM, s=15., marker='o', c='g', linewidths=[0.,], zorder=0)
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



def puls_heatmap(ax, puls, idL, db, pulseN=False, inc=12):
  if inc == 12:    
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
  else:
    ra = np.array([ 497,  497,  582,  582,  497,  417,  417,  497,  582,  668,  668,
        668,  582,  497,  417,  331,  331,  331,  417,  497,  582,  668,
        748,  748,  748,  748,  668,  582,  497,  417,  331,  251,  251,
        251,  251,  331,  417,  497,  582,  668,  748,  834,  834,  834,
        834,  834,  748,  662,  582,  497,  417,  337,  251,  165,  165,
        165,  165,  165,  251,  331,  417,  497,  582,  668,  748,  834,
        914,  914,  914,  914,  914,  914,  828,  748,  662,  582,  497,
        417,  337,  251,  171,   85,   85,   85,   85,   85,   85,  165,
        251,  331,  417,  497,  582,  668,  748,  834,  914, 1000, 1000,
       1000, 1000, 1000, 1000,  994,  914,  828,  748,  662,  582,  497,
        417,  337,  251,  171,   85,    5,    0,    0,    0,    0,    0,
          0,   85,  165,  251,  331,  417])
    dec = np.array([ 499,  583,  541,  458,  416,  458,  541,  666,  624,  583,  499,
        416,  375,  333,  375,  416,  499,  583,  624,  750,  708,  666,
        624,  541,  458,  375,  333,  291,  249,  291,  333,  375,  458,
        541,  624,  666,  708,  833,  791,  750,  708,  666,  583,  499,
        416,  333,  291,  249,  208,  166,  208,  249,  291,  333,  416,
        499,  583,  666,  708,  750,  791,  916,  874,  833,  791,  750,
        708,  624,  541,  458,  375,  291,  249,  208,  166,  125,   83,
        125,  166,  208,  249,  291,  375,  458,  541,  624,  708,  750,
        791,  833,  874, 1000,  957,  916,  874,  833,  791,  750,  666,
        583,  499,  416,  333,  249,  208,  166,  125,   83,   42,    0,
         42,   83,  125,  166,  208,  249,  333,  416,  499,  583,  666,
        750,  791,  833,  874,  916,  957])
  
  n_beams = ra.size
  dDM = 0.2
  dm_l = float(puls.DM - dDM/2. - 0.001)
  dm_h = float(puls.DM + dDM/2. + 0.001)
  t_l = float(puls.Time - 2. * puls.Duration)
  t_h = float(puls.Time + 2. * puls.Duration)
  sap = int(puls.SAP)
  select = '(SAP == sap) and ((DM > dm_l) and (DM < dm_h)) and ((Time >= t_l) and (Time <= t_h))'
  try:
    events = pd.read_hdf(db, 'events', where=select)
  except ValueError:
    events = pd.read_hdf(db, 'events')
    events = events[(events.SAP == sap) & ((events.DM > dm_l) & (events.DM < dm_h)) & ((events.Time >= t_l) & (events.Time <= t_h))]
  SNR = events.groupby('BEAM').Sigma.sum()
  ind = pd.Series(np.zeros(n_beams))
    
  ind.index += inc + 1
  SNR = SNR.reindex_like(ind)

  plot = ax.scatter(ra,dec,s=200,edgecolor='none',c=SNR,cmap='hot_r')
  
  #bar = plt.colorbar(plot, ax=ax)
  #bar.set_label('Cumulative SNR')
  
  ax.set_xlabel('RA')
  ax.set_ylabel('DEC')
  if pulseN: 
    ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN))
    ax.annotate('DM: {:.2f}$\pm$0.2, Time: {:.2f}$\pm${:.2f}'.format(puls.DM,puls.Time,puls.Duration), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  else:
    ax.annotate('DM: {:.2f} - {:.2f}, Time: {:.2f} - {:.2f}'.format(dm_l,dm_h,t_l,t_h), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  
  beam = int(puls.BEAM)
  if beam > inc: ax.scatter(ra[beam-(inc+1)],dec[beam-(inc+1)],s=300,linewidths=[0.,],marker='*',c='w')
  [ax.annotate(str(i+(inc+1)),(ra[i],dec[i]), horizontalalignment='center', verticalalignment='center', color='m') for i in range(0,n_beams)]

  ax.set_xlim(-100,1100)
  ax.set_ylim(-100,1100)
  ax.tick_params(axis='both', labelleft='off', labelright='off', labeltop='on', labelbottom='off', right='off', left='off', bottom='off', top='off')
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

  ax.axis('off')
  return



def puls_dynSpec(ax1, ax2, puls, idL, inc=12):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam == inc: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  if inc == 0: conf = 'confirmations'
  else: conf = ''
  filename = '{folder}/{conf}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  maskfn = '{folder}/{conf}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}_rfifind.mask'.format(folder=RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(filename): return -1
  if os.path.isfile(maskfn): mask = True  
  else: mask = False
  
  duration = 10
  
  filetype, header = Utilities.read_header(filename)
  MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
  v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
  
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits.PsrfitsFile(filename), puls.Time_org*(1+v)-puls.Duration*duration/2, puls.Duration*(duration+1), dm=puls.DM, downsamp=int(puls.Downfact), nsub=16, maskfn=maskfn, mask=mask)
  waterfaller.plot_waterfall(ds, start, puls.Duration*duration, ax_im=ax1, interactive=False)
  ax1.scatter((start + puls.Duration)/2.,F_MIN,marker='^',s=25,c='r',lw=0.)
  
  DM_delay = presto.psr_utils.delay_from_DM(puls.DM, F_MIN) - presto.psr_utils.delay_from_DM(puls.DM, F_MAX)
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits.PsrfitsFile(filename), puls.Time_org*(1+v)-0.1, puls.Duration+DM_delay+0.2, downsamp=int(puls.Downfact), nsub=32, maskfn=maskfn, mask=mask)
  waterfaller.plot_waterfall(ds, start, puls.Duration+DM_delay+0.1, ax_im=ax2, interactive=False, sweep_dms=[puls.DM])

  #ax2.set_xlabel('Time (s)')
  #ax2.set_ylabel('Frequency (MHz)')  
  #ax1.set_xlabel('Time (s)')
  #ax1.set_ylabel('Frequency (MHz)')
      
  return 0



def load_ts(puls, idL, filename):
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

  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  out_dir = TMP_FOLDER.format(idL) + '/timeseries/'
  if not os.path.isdir(os.path.join(TMP_FOLDER.format(idL), 'timeseries')): os.makedirs(os.path.join(TMP_FOLDER.format(idL), 'timeseries'))
  FNULL = open(os.devnull, 'w')
  for j,DM in enumerate(DM_range):
    if not os.path.isfile(filename.format(DM)):
      error = subprocess.call(['sh', '/projects/0/lotaas2/software/LSP/spdspsr_pl.sh', idL, str(sap), str(beam), '{:.2f}'.format(DM), out_dir], stdout=FNULL, stderr=FNULL)
    
    try:
      ts = np.memmap(filename.format(DM), dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
    except IOError: ts = np.zeros(nBins) + np.nan

    data[j] = ts
    
  params = {'bins_out': nBins*scrunch_fact, 'bin_start': bin_start, 'scrunch_fact': scrunch_fact, 'duration': duration, 'DM_min': DM_range[0], 'DM_max': DM_range[-1], 'k': k}
  
  return data, params



def puls_dedispersed(ax, puls, idL, pulseN=False, inc=12):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam == inc: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  if inc == 0: conf = 'confirmations'
  else: conf = ''
  raw_dir = '{folder}/{conf}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(raw_dir): return -1
  
  filename = TMP_FOLDER.format(idL) + '/timeseries/{}/manual_fold_DM{{0:.2f}}.dat'.format(idL)
  data, params = load_ts(puls, idL, filename)

  #Image plot
  ax.imshow(data,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=[0,params['bins_out'],params['DM_min'], params['DM_max']])
  ax.set_ylim((params['DM_min'], params['DM_max']))
  ax.set_ylabel('DM (pc/cc)')
  if pulseN: ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN), y=1.08)

  #Plot contours
  x = np.arange(params['bins_out'])
  y = (x - params['bins_out'] / 2.) / params['k'] + puls.DM
  ax.plot(x, y,'r--')
  ax.axvline((params['bins_out'] - params['scrunch_fact']) / 2., color='r', ls='--')

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
      ts = np.memmap(filename.format(puls.DM), dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      ts = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
    except IOError: x = ts = None
    return x, ts
  x, ts = inset(filename, params, puls)
  
  ax3 = inset_axes(ax, width="30%", height="30%", loc=1)
  ax3.plot(x, ts, 'k')
  ax3.set_xticks([])
  ax3.set_yticks([])

  #Time axis
  ax.set_xlim((0,params['bins_out']))
  labels = ax.get_xticks().tolist()
  if puls.DM < 40.52: down = 1
  elif DM_peak < 141.77: down = 2
  else: down = 4
  new_labels = ["%.3f" % (n + params['bin_start'] * RES * down) for n in labels]
  ax.set_xticklabels(new_labels)
  ax.set_xlabel('Time (s)')

  #Sample axis
  ax2 = ax.twiny()
  ax2.set_xlim((0,params['bins_out']))
  labels = ax.get_xticks().tolist()
  new_labels = [int((n + params['bin_start']) % 100000) for n in labels]
  ax2.set_xticklabels(new_labels)
  ax2.set_xlabel('Sample - {}00000'.format(params['bin_start'] / 100000))
  

  return 0
