#!/usr/bin/env python

import os
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import presto
import pandas as pd
import pyfits
import pyximport; pyximport.install()

import C_Funct
import Utilities

#Install Cfunct if not present
if not os.path.isfile('repeating_FRB_Cfunct'):
  from distutils.core import setup
  from distutils.extension import Extension
  from Cython.Distutils import build_ext
  ext_modules = [Extension("src.C_Funct", ["src/C_Funct.pyx"],include_dirs=[np.get_include()])]
  setup(name = 'LSPs',cmdclass = {'build_ext': build_ext},ext_modules = ext_modules)
  #NOT WORKING: find a way to give python setup.py build_ext --inplace
import repeating_FRB_Cfunct as C_Funct


singlepulse_filename = '/{}'
fits_filename = '/{}'
plot_folder = '/{}'
N_beams = 7

DM_SPACING = 
dT_GROUP = 
dDM_GROUP = 
dt = 0.01 #s #Minimum repetition tim


inf = read_header()
F_MIN = inf['OBSFREQ'] - inf['CHAN_BW'] * inf['NCHAN'] / 2.
F_MAX = inf['OBSFREQ'] + inf['CHAN_BW'] * inf['NCHAN'] / 2.
RES = inf['TBIN']

def read_header():
  for beam in N_beams:
    if os.path.isfile(fits_filename.format(beam)):
      fits = pyfits.open(fits_filename.format(beam),memmap=True)
      header = fits['SUBINT'].header + fits['PRIMARY'].header
      fits.close()
      return header
  print 'No fits file found - exiting...'
  exit()


def load_candidates():
  def read_candidates():
    def load_pulses():    
      def load_beam_pulses(beam):
        def read_singlepulse(beam):
          filename = singlepulse_filename.format(beam)
          events = pd.read_csv(singlepulse_file, delim_whitespace=True)
          events.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']  #TO MODIFY
          events['Duration'] = events.Sampling * events.Downfact
          events = events.ix[:,['DM','Sigma','Time','Sample','Downfact','Sampling']]
          events.index.name = 'idx'
          events['BEAM'] = beam
          events['Pulse'] = 0

        events = read_singlepulse(beam)
          
        def straight_time(Time, DM):
          k = 4149. #s-1
          delay = k * (F_MIN**-2 - F_MAX**-2)
          Time += np.float32( delay * DM / 2)
          return Time
          
        events['Time'] = straight_time(events.Time,events.DM)
          
        def group_dm(events):
          events.sort(['DM','Time'],inplace=True)
          C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Duration.values, events.Pulse.values, DM_SPACING, dT_GROUP, dDM_GROUP)
          events = events[events.Pulse >= 0]
          events.Pulse = events.Pulse * 10 + events.BEAM
          gb = events.groupby('Pulse',sort=False)
          pulses = events.loc[gb.Sigma.idxmax()]  
          pulses.index = pulses.Pulse
          pulses.index.name = 'idx'
          pulses['Candidate'] = -2
          pulses['N_events'] = gb.DM.count()
          return pulses

        pulses = group_dm(events)
        pulses = pulses[pulses.N_events > 2]
        return pulses

      pool = mp.Pool()
      results = pool.map_async(load_beam_pulses, range(N_beams))  #TO CHECK
      pool.close()
      pool.join()
      return pd.concat([p.get() for p in results])
    
    pulses = load_pulses()
      
    def group_beam():
      pulses.sort(['Time','Sigma'],inplace=True)
      C_Funct.Get_Candidate(events.Beam.values, events.Sigma.values, events.Time.values, events.Candidate.values, events.Duration.values)
      if pulses[pulses.Candidate < -1].shape[0] > 0: print "ATTENTION! - Bug in Cfunct!"
      pulses = pulses[pulses.Candidate >= 0]
      gb = pulses.groupby('Candidate',sort=False)
      candidates = pulses.loc[gb.Sigma.idxmax()]  
      candidates.index = candidates.Candidate
      candidates.index.name = 'idx'
      candidates['N_pulses'] = gb.DM.count()

      def rfi_filters(candidates,gb):
        #Minimum sigma ratio for pulses appearing in more than 3 beams
        Sigma_ratio = gb.Sigma.max() / gb.Sigma.min()
        candidates = candidates[(candidates.N_pulses <= 3) | (Sigma_ratio > 1.2)]
        return candidates
      
      candidates = rfi_filters(candidates,gb)
      
      return candidates
    
    candidates = group_beam(pulses)
    return candidates
    
  candidates = read_candidates()
  candidates.sort('Sigma',inplace=True)
  
  return candidates.head(20)
  
  
  

def plot_spectrum(cands):
  def plot(cands):
    for puls in cands.iterrows():
      beam = puls[1]['BEAM']

      duration = np.int(np.round(puls.Duration/RES))
      spectra_border = 20
      offset = duration*spectra_border
      sample = puls.Sample * puls.Downfact
      
      freq = np.linspace(F_MIN,F_MAX,2592)
      time = (4149 * puls.DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
      
      def bary(sample):
        MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
        v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
        sample += np.round(sample*v).astype(int)
        return sample
      sample = bary(sample)
      
      def load_spectra(filename,mask_file):
        def dedispersion(spectrum):  
          for i in range(spectrum.shape[1]):
            spectrum[:,i] = np.roll(spectrum[:,i], time[i])
          spectrum = spectrum[:2*offset+duration]
          return spectrum

        spectrum = Utilities.read_fits(filename,puls.DM.copy(),sample.copy(),duration,offset,RFI_reduct=True,mask=mask_file)
        dedisp_spectrum = dedispersion(spectrum.copy())

        spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0]/duration,duration,spectrum.shape[1]]),axis=1)
        spectrum = np.mean(np.reshape(spectrum,[spectrum.shape[0],spectrum.shape[1]/32,32]),axis=2)
        dedisp_spectrum = np.mean(np.reshape(dedisp_spectrum,[dedisp_spectrum.shape[0]/duration,duration,dedisp_spectrum.shape[1]]),axis=1)
        dedisp_spectrum = np.mean(np.reshape(dedisp_spectrum,[dedisp_spectrum.shape[0],dedisp_spectrum.shape[1]/32,32]),axis=2)
        
        return spectrum, dedisp_spectrum
      
      spectrum, dedisp_spectrum = load_spectra(fits_file.format(beam),mask_file.format(beam))
      
      plt.clf()
      #fig = plt.figure()

      ax1 = plt.subplot2grid((5,3),(0,0), rowspan=4, colspan=2) #raw spectrum
      ax2 = plt.subplot2grid((5,3),(1,2), rowspan=3, colspan=1) #dedisp spectrum
      ax3 = plt.subplot2grid((5,3),(0,2), rowspan=1, colspan=1) #dedisp timeseries
      ax4 = plt.subplot2grid((5,3),(4,0), rowspan=5, colspan=1) #meta data 
      
      def plot_raw_spectrum(ax,spectrum):
        extent = [(sample-offset)*RES,(sample-offset+spectrum.shape[0])*RES,F_MIN,F_MAX]
        ax.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
        ax.plot(-time, freq, 'r--')
        ax.plot(-time+2*offset, freq, 'r--')
        ax.axis(extent)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Frequency (MHz)')
        
      def plot_dedisp_spectrum(ax,spectrum):
        extent = [(sample-offset)*RES,(sample+duration+offset)*RES,F_MIN,F_MAX]
        ax.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
        ax.scatter((sample+duration/2)*RES,F_MIN+1,marker='^',s=1000,c='r')
        ax.axis(extent)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Frequency (MHz)') 
        
      def plot_timeseries(ax,spectrum):
        timeseries = np.sum(spectrum, axis=1)
        ax.plot(timeseries, 'k-')
        ax.xaxis.set_major_formatter(plt.NullFormatter())
      
      plot_raw_spectrum(ax1,spectrum)
      plot_dedisp_spectrum(ax2,dedisp_spectrum)
      plot_timeseries(ax3,dedisp_spectrum)
      plot_metadata(ax4)
      
      plt.savefig(plot_folder.format(beam),format='png',bbox_inches='tight',dpi=200)
      
  
  def beams_parallel(cands):
    CPUs = mp.cpu_count()
    cands_per_CPU = int(np.ceil(cands.shape[0]/CPUs))
    pool = mp.Pool()
    [pool.apply_async(plot, args=cands.iloc[i:(i+1)*cands_per_CPU]) for i in range(CPUs)]
    pool.close()
    pool.join()
    return

  beams_parallel(cands)
  

  
def main():
  candidates = load_candidates()  
  plot_spectrum(candidates)


if __name__ == "__main__":
  main()
  
 
  
  
  