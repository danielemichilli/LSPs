#############################
# LOTAAS Single Pulse plots
# Written by Daniele Michilli
#############################

from glob import glob
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
from waterfaller import waterfaller
import psrfits
import re

from Parameters import *
import RFIexcision
import Paths as PATH

mpl.rc('font', size=4)
mpl.rc('lines', linewidth=1)


def output(idL, pulses, meta_data, candidates, db, inc=12):
  pulses.sort_values('Sigma',ascending=False,inplace=True)
  candidates.sort_values('Sigma',ascending=False,inplace=True)
  
  plt.close('all')
  fig = plt.figure(figsize=(8,4))

  out_dir = os.path.join(PATH.TMP_FOLDER, 'timeseries')
  if os.path.isdir(out_dir): shutil.rmtree(out_dir)
  for idx_c, cand in candidates.iterrows():
    sap = int(cand.SAP)
    beam = int(cand.BEAM)
    events = pd.read_hdf(PATH.DB,'events',where=['(SAP == sap) & (BEAM == beam)'])
    store = os.path.join(PATH.WRK_FOLDER, 'sp/candidates/{}.pdf'.format(cand.id))
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

  gs = gridspec.GridSpec(3, 8, wspace=.8, hspace=.5)
  ax1 = plt.subplot(gs.new_subplotspec((0,0), 2, 2), rasterized = True)
  ax2 = plt.subplot(gs.new_subplotspec((0,2), 2, 6), rasterized = True)
  ax3 = plt.subplot(gs.new_subplotspec((2,0), 1, 2), rasterized = True)
  ax4 = plt.subplot(gs.new_subplotspec((2,2), 1, 2), rasterized = True)
  ax5 = plt.subplot(gs.new_subplotspec((2,4), 1, 2), rasterized = True)
  ax6 = plt.subplot(gs.new_subplotspec((2,6), 1, 2), rasterized = True)
  
  meta_data_plot(ax1,meta_data[(meta_data.SAP==sap)&(meta_data.BEAM==beam)],pulses,cand)
  scatter_beam(ax2, pulses, pulses_beam, cand)

  if not pulses_beam.empty: 
    scatter_SNR(ax4,pulses,events[events.Pulse.isin(pulses_beam.index)],cand)
    scatter_SNR(ax5,pulses,events[events.Pulse.isin(pulses_beam.index)],cand, zoom=True)
  try:
    hist_DM(ax3,pulses_beam,cand)
    hist_SNR(ax6,pulses_beam,cand)
  except ValueError: pass
  
  pdf.savefig(bbox_inches='tight', dpi=200)
  return



def puls_plot(pdf, puls, events, idL, db, i, inc=12):
  gs = gridspec.GridSpec(2, 6, wspace=0.5, hspace=0.2)
  ax1 = plt.subplot(gs.new_subplotspec((0,0), 1, 1), rasterized = True)
  ax2 = plt.subplot(gs.new_subplotspec((0,1), 1, 1), rasterized = True)
  ax3 = plt.subplot(gs.new_subplotspec((1,0), 1, 2), rasterized = True)
  ax4 = plt.subplot(gs.new_subplotspec((0,2), 1, 1), rasterized = True)
  ax5 = plt.subplot(gs.new_subplotspec((0,3), 1, 1), rasterized = True)
  ax6 = plt.subplot(gs.new_subplotspec((1,2), 1, 1), rasterized = True)
  ax7 = plt.subplot(gs.new_subplotspec((1,3), 1, 1), sharey=ax6, rasterized = True)
  ax8 = plt.subplot(gs.new_subplotspec((0,4), 2, 2), rasterized = True)

  ev = events[events.Pulse == puls.name]

  puls_meta_data(ax1, puls, ev.Pulse.iloc[0], i)
  puls_DM_Time(ax2, ev, events, puls)
  puls_dedispersed_plotted = puls_dedispersed(ax3, puls, idL, inc=inc, prof_ax=ax4)
  if puls_dedispersed_plotted == -1: 
    plot_not_valid(ax3)
    puls_dynSpec_plotted = puls_dynSpec(ax6, ax7, puls, idL, inc=inc, ax_ts=ax4)
  else:
    puls_dynSpec_plotted = puls_dynSpec(ax6, ax7, puls, idL, inc=inc, ax_ts=None)
  if puls_dynSpec_plotted == -1:
    plot_not_valid(ax6)
    plot_not_valid(ax7)
    if puls_dedispersed_plotted == -1: 
      plot_not_valid(ax4)
  puls_SNR_DM(ax5, ev)
  puls_heatmap(ax8, puls, idL, db, inc=inc)

  pdf.savefig(bbox_inches='tight', dpi=200)
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
    new_val = np.clip(values, 5, 20)
    m = 31.85
    q = -137.025
    return new_val * m + q
  sig = circle_size(pulses_beam.Sigma)

  ax.scatter(pulses_beam.Time, pulses_beam.DM, c='k', s=sig, edgecolor='w', lw=.2, zorder=2)

  ax.axhline(cand.DM, c='r', ls='-', zorder=1)
  
  ax.axhline(DM_STEP1,c='k',ls='--',lw=.1, zorder=1)
  ax.axhline(DM_STEP2,c='k',ls='--',lw=.1, zorder=1)
  ax.set_yscale('log')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('DM (pc/cc)')
  ax.axis([0,3600,DM_MIN,550])
  
  sig = circle_size(pulses.Sigma)
  ax.scatter(pulses.Time, pulses.DM, s=sig, linewidths=[0.,], marker='*', c='w', zorder=3)
  
  top = pulses.iloc[:10]
  for i in range(top.shape[0]):
    ax.annotate(str(i),xy=(top.Time.iloc[i],top.DM.iloc[i]),horizontalalignment='center',\
      verticalalignment='center',color='dodgerblue',size=5,fontweight='bold',zorder=4)
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
    ax.annotate('Width (ms): {0:.0f}'.format(pulses.Duration.iloc[0]*1000), xy=(0,2))
    ax.annotate('N. events: {}'.format(pulses.N_events.iloc[0]), xy=(0,1))
  else:
    ax.annotate('dDM (pc/cm2): {0:.2f}'.format(pulses.DM.max()-pulses.DM.min()), xy=(0,4))
    ax.annotate('Period (s): {0:.3f}'.format(cand.Period), xy=(0,3))
    ax.annotate('Period err. (s): {0:.3f}'.format(cand.Period_err), xy=(0,2))
    ax.annotate('N. pulses: {0:.0f}'.format(cand.N_pulses), xy=(0,1))
  ax.annotate('L-SpS {}'.format(meta_data.version.iloc[0]), xy=(0,0))
  ax.axis('off')
  return



def hist_DM(ax,pulses,cand):
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)
  ax.hist(pulses.DM.tolist(),bins=int(550-DM_MIN),histtype='stepfilled',color=u'k',range=(DM_MIN,550), zorder=2)
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('Counts')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(DM_STEP1,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(DM_STEP2,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def hist_SNR(ax,pulses,cand):
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)
  ax.hist(pulses.DM.tolist(),bins=int(550-DM_MIN),histtype='stepfilled',color=u'k',weights=pulses.Sigma.tolist(),range=(DM_MIN,550), zorder=2)
  ax.set_xscale('log')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('Cumulative S/N')
  ax.set_xlim(DM_MIN,550)
  ax.axvline(DM_STEP1,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(DM_STEP2,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def scatter_SNR(ax, pulses, pulses_beam, cand, zoom=False):
  ax.axvline(cand.DM, c='r', ls='-', linewidth=.2, zorder=3)
  ax.scatter(pulses_beam.DM,pulses_beam.Sigma,c='k',s=5,linewidths=[0.,], zorder=2)
  if not zoom: ax.set_xscale('log')
  ax.set_ylabel('S/N')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  if zoom: ax.set_xlim((cand.DM - cand.dDM * 2, cand.DM + cand.dDM * 2))
  else: ax.set_xlim((3., 550.))
  ax.set_ylim((6.5, pulses_beam.Sigma.max()+1))
  ax.axvline(DM_STEP1,c='k',ls='--',lw=.1, zorder=1)
  ax.axvline(DM_STEP2,c='k',ls='--',lw=.1, zorder=1)
  ax.tick_params(which='both',direction='out')
  ax.yaxis.set_ticks_position('left')
  return



def puls_DM_Time(ax, event, all_events, puls):
  def circle_size(values):
    new_val = np.clip(values,5,20)
    m = 31.85
    q = -137.025
    return (new_val * m + q) / 2. + 5.
  sig = circle_size(event.Sigma)
  #sig = event.Sigma - event.Sigma.min()
  #sig = sig * (480. / sig.max()) + 20.
  
  x_ev = (event.Time - puls.Time)
  ax.scatter(x_ev, event.DM, s=sig, marker='o', linewidths=[.5,], edgecolor='r', facecolors='none', zorder=0)
  #ax.scatter(event.Time, event.DM, s=20., marker='o', c='k', linewidths=[0.,], zorder=1)
  #ax.errorbar(puls.Time, puls.DM, xerr=puls.Duration/2, yerr=puls.dDM/2, lw=.2, fmt='none', ecolor='r', zorder=2)
  ax.scatter(0, puls.DM, marker='x', s=30., color='dodgerblue', zorder=2)
  ax.set_xlabel('$\Delta$Time (s)')
  ax.set_ylabel('DM (pc cm$^{-3}$)')
  
  DM_min = np.clip(puls.DM - puls.dDM * 3., 3., 550.)
  DM_max = np.clip(puls.DM + puls.dDM * 3., 3., 550.)
  k = 4148.808 #s-1
  delay = k * (F_MIN**-2 - F_MAX**-2)
  dDM = DM_max - DM_min
  dt = delay * dDM * 3.
  Time_min = puls.Time - dt
  Time_max = puls.Time + dt
  events = all_events[(all_events.DM > DM_min) & (all_events.DM < DM_max) & (all_events.Time > Time_min) & (all_events.Time < Time_max)]
  x_ev = (events.Time - puls.Time)
  ax.scatter(x_ev, events.DM, s=15., marker='o', c='k', linewidths=[0.,], zorder=1)
  return



def puls_SNR_DM(ax, event):
  ax.plot(event.DM, event.Sigma, 'k')
  ax.set_xlabel('DM (pc cm$^{-3}$)')
  ax.set_ylabel('S/N')
  
  ax2 = ax.twinx()
  ax2.plot(event.DM, event.Duration*1000, 'r')
  ax2.set_ylabel('Width (ms)', color='r')
  for tl in ax2.get_yticklabels():
    tl.set_color('r')

  ax.set_zorder(ax2.get_zorder()+1)
  ax.patch.set_visible(False)
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

  '''
  dm_l = float(puls.DM - 2/2. - 0.001)
  dm_h = float(puls.DM + 2/2. + 0.001)
  t_l = float(puls.Time - 100. * puls.Duration)
  t_h = float(puls.Time + 100. * puls.Duration)
  '''

  sap = int(puls.SAP)
  select = '(SAP == sap) and ((DM > dm_l) and (DM < dm_h)) and ((Time >= t_l) and (Time <= t_h))'
  try:
    events = pd.read_hdf(db, 'events', where=select)
  except (TypeError, ValueError) as e:
    events = pd.read_hdf(db, 'events')
    events = events[(events.SAP == sap) & ((events.DM > dm_l) & (events.DM < dm_h)) & ((events.Time >= t_l) & (events.Time <= t_h))]
  SNR = events.groupby('BEAM').Sigma.max() #.sum()

  ind = pd.Series(np.zeros(n_beams))
    
  ind.index += inc + 1
  SNR = SNR.reindex_like(ind)

  SNR -= SNR.min()
  SNR /= SNR.max()

  plot = ax.scatter(ra,dec,s=200,edgecolor='none',c=SNR,cmap='hot_r')
  
  bar = plt.colorbar(plot, ax=ax, orientation='horizontal', pad=0.03)# ,ticks=[])  
  bar.set_label('Normalized S/N')
  
  ax.set_xlabel('RA')
  ax.xaxis.set_label_position('top')
  ax.set_ylabel('DEC')
  if pulseN: 
    ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN))
    ax.annotate('DM: {:.2f}$\pm$0.2, Time: {:.2f}$\pm${:.2f}'.format(puls.DM,puls.Time,puls.Duration), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  else:
    ax.annotate('DM: {:.2f} - {:.2f}, Time: {:.2f} - {:.2f}'.format(dm_l,dm_h,t_l,t_h), xy=(-80,1080), fontsize='large',horizontalalignment='left',verticalalignment='top')
  
  beam = int(puls.BEAM)
  if beam > inc: ax.scatter(ra[beam-(inc+1)],dec[beam-(inc+1)],s=300,linewidths=[0.,],marker='*',c='w')
  [ax.annotate(str(i+(inc+1)),(ra[i],dec[i]), horizontalalignment='center', verticalalignment='center', color='dodgerblue') for i in range(0,n_beams)]

  ax.set_xlim(-100,1100)
  ax.set_ylim(-100,1100)
  ax.tick_params(axis='both', labelbottom='off', labelleft='off', top='off', bottom='off', left='off', right='off')
  return



def puls_meta_data(ax, puls, idx, i):
  ax.axis([0,10,0,10])
  ax.annotate("Pulse n. {}".format(i), xy=(0,9)) 
  ax.annotate('Pulse code: {}'.format(idx), xy=(0,8))
  ax.annotate('DM (pc/cm2): {0:.2f}'.format(puls.DM), xy=(0,7))
  ax.annotate('dDM (pc/cm2): {0:.2f}'.format(puls.dDM), xy=(0,6))
  ax.annotate('Time (s): {0:.2f}'.format(puls.Time_org), xy=(0,5))
  ax.annotate('Sigma: {0:.1f}'.format(puls.Sigma), xy=(0,4))
  ax.annotate('Width (ms): {0:.0f}'.format(puls.Duration*1000), xy=(0,3))
  ax.annotate('N events: {0:.0f}'.format(puls.N_events), xy=(0,2))
  ax.annotate('Sample: {}'.format(int(np.round(puls.Time_org / RES))), xy=(0,1))
  if puls.DM < DM_STEP1: d = 1
  elif puls.DM < DM_STEP2: d = 2
  else: d = 4
  ax.annotate('Downsampling: {}'.format(puls.Downfact * d), xy=(0,0))
  ax.axis('off')
  return



def puls_dynSpec(ax1, ax2, puls, idL, inc=12, ax_ts=None):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam == inc: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  if inc == 0: conf = 'confirmations'
  else: conf = ''
  filename = '{folder}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=PATH.RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  print filename
  maskfn = '{folder}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}_rfifind.mask'.format(folder=PATH.RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  print maskfn
  if not os.path.isfile(filename): return -1
  if os.path.isfile(maskfn): mask = True  
  else: mask = False
  
  duration = 20

  try: df = int(puls.Downfact)
  except AttributeError:
    if puls.DM < DM_STEP1: df = int(round(puls.Duration / 0.0004915))
    elif puls.DM < DM_STEP2: df = int(round(puls.Duration / 0.0004915 / 2))
    else: df = int(round(puls.Duration / 0.0004915 / 4))
  
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits.PsrfitsFile(filename), puls.Time_org-puls.Duration*duration/2, puls.Duration*(duration+1), nsub=16, dm=puls.DM, width_bins=df, maskfn=maskfn, mask=mask, scaleindep=False, bandpass_corr=True, zerodm=True)
  waterfaller.plot_waterfall(ds, start, puls.Duration*(duration+1), ax_im=ax1, interactive=False, puls_t=-puls.Duration*duration/2, ax_ts=ax_ts)
  ax1.scatter(0.,F_MIN+1,marker='^',s=50,c='r',lw=0.)

  DM_delay = presto.psr_utils.delay_from_DM(puls.DM, F_MIN) - presto.psr_utils.delay_from_DM(puls.DM, F_MAX)
  ds, nbinsextra, nbins, start = waterfaller.waterfall(psrfits.PsrfitsFile(filename), puls.Time_org-0.1, puls.Duration+DM_delay+0.2, dm=0., nsub=32*3, downsamp=df*3, maskfn=maskfn, mask=mask, scaleindep=True, zerodm=True, bandpass_corr=True)
  waterfaller.plot_waterfall(ds, start, puls.Duration+DM_delay+0.2, ax_im=ax2, interactive=False, sweep_dms=[puls.DM], puls_t=-0.1)

  return 0



def load_ts(puls, idL, filename):
  out_dir = os.path.join(PATH.TMP_FOLDER, 'timeseries')
  try: shutil.rmtree(out_dir)
  except OSError: pass
  os.makedirs(out_dir)

  nDMs = 21
  DM = puls['DM']
  dDM = puls['dDM'] * 4
  lowDM = DM - dDM / 2.
  stepDM = dDM / (nDMs - 1)
  
  mask = os.path.splitext(filename)[0] + "_rfifind.mask"
  error = subprocess.call(['prepsubband', '-dmprec', '4', '-numdms', str(nDMs), '-dmstep', str(stepDM), \
    '-nsub', '288', '-lodm', str(lowDM), '-mask', mask, '-runavg', '-noscales', '-noweights',\
    '-nooffsets', '-o', 'diagnostic_plot', filename], cwd=out_dir)
      
  nProfBins = 3
  duration = int(np.round(puls['Duration'] / RES))
  k = 4148.808 * (F_MIN**-2 - F_MAX**-2) / RES
  nPlotBins = int(np.ceil(dDM * k / duration))
  nBins = nPlotBins * nProfBins
  data = np.zeros((nDMs,nBins))
  peak = int(np.round(puls['Time_org'] / RES))
  scrunch_fact = int(np.round(duration / float(nProfBins)))
  if scrunch_fact < 1: scrunch_fact = 1 
  bin_start = peak - (nBins * scrunch_fact - 2) / 2

  def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
    
  ts_list = natural_sort(glob(os.path.join(out_dir, 'diagnostic_plot*.dat')))
  for i,ts_name in enumerate(ts_list):
    try:
      ts = np.memmap(ts_name, dtype=np.float32, mode='r', offset=bin_start*4, shape=(nBins*scrunch_fact,))
      data[i] = np.mean(np.reshape(ts, (nBins, scrunch_fact)), axis=1)
    except IOError: data[i] = np.zeros(nBins) + np.nan  
  
  idx = len(ts_list) / 2
  prof = np.fromfile(ts_list[idx], dtype=np.float32)
  return data, nBins*scrunch_fact*RES, prof



def puls_dedispersed(ax, puls, idL, pulseN=False, inc=12, prof_ax=False):
  sap = int(puls.SAP)
  beam = int(puls.BEAM)
  if beam == inc: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  if inc == 0: conf = 'confirmations'
  else: conf = ''
  raw_dir = '{folder}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=PATH.RAW_FOLDER,conf=conf,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(raw_dir): return -1

  data, plot_duration, prof = load_ts(puls, idL, raw_dir)

  #Image plot
  ax.imshow(data, cmap='Greys', origin="lower", aspect='auto', interpolation='nearest',\
    extent=[-plot_duration / 2., plot_duration / 2., puls.DM - puls.dDM / 2., puls.DM + puls.dDM / 2.])
  ax.set_ylim((puls.DM - puls.dDM / 2., puls.DM + puls.dDM / 2.))
  ax.set_ylabel('DM (pc cm$^{-3}$)')
  if pulseN: ax.set_title('{obs} SAP{sap} BEAM{beam} - Candidate {cand} Pulse {puls}'.format(obs=idL,sap=puls.SAP,beam=puls.BEAM,cand=puls.Candidate,puls=pulseN), y=1.08)

  #Plot contours
  k = 4148.808 * (F_MAX**-2 - F_MIN**-2)
  x = np.array([-plot_duration/2., plot_duration/2.])
  y = x / k + puls.DM
  ax.plot(x, y,'r', linewidth=.2)
  ax.axvline(0, color='r', linewidth=.2)

  def inset(prof):
    idx = int(np.round(puls.Time_org / RES))
    prof_bins = 4
    scrunch = int(puls.Downfact / prof_bins)
    if scrunch < 1: scrunch = 1
    bins = 101
    prof = prof[idx - bins / 2 * scrunch : idx + bins / 2 * scrunch]
    prof = prof.reshape([bins, scrunch]).sum(axis=1)
    x = np.linspace(-(bins/2), bins/2, bins) * RES * 1e3 * prof_bins
    return x, prof

  ts = inset(prof)
  
  if not prof_ax: prof_ax = inset_axes(ax, width="30%", height="30%", loc=1)
  prof_ax.plot(ts, 'k')
  prof_ax.set_xlim((0,ts.size))
  prof_ax.set_yticks([])
  prof_ax.set_xlabel('Time (ms)')

  #Time axis
  ax.set_xlim((-plot_duration/2.,plot_duration/2.))
  if puls.DM < DM_STEP1: down = 1
  elif puls.DM < DM_STEP2: down = 2
  else: down = 4
  ax.set_xlabel('$\Delta$Time (s)') #.format(params['bin_start'] * RES * down * 1000 * params['scrunch_fact']))

  return 0
