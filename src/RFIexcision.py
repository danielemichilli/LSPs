########################################
#
# Radio Frequency Interferences excision
#
# Written by Daniele Michilli
#
########################################

import numpy as np
import pandas as pd
import logging
from scipy import special
from scipy import stats
import os
import subprocess
try: 
  import pyfits
  import presto
except ImportError: pass

import Utilities
import C_Funct
import Paths
from Parameters import *


def sift_pulses(pulses, events, idL, sap, beam):
  arff_basename = '{}/thresholds_{}_{}'.format(Paths.TMP_FOLDER.format(idL), sap, beam)
  filters(pulses, events, arff_basename+'.arff')
  ML_predict = os.path.join(Paths.TMP_FOLDER.format(idL), 'ML_predict.txt')  
  pulses = select_real_pulses(pulses,arff_basename, ML_predict)
  return pulses
  
  

def select_real_pulses(pulses,basename, out_name):
  subprocess.call(['java', '-jar', Paths.CLASSIFIER, '-v', '-m{}'.format(Paths.MODEL_FILE), '-p{}'.format(basename+'.arff'), '-o{}'.format(basename+'.positive'), '-a1'])
  os.remove(basename+'.arff')
  pulses_list = np.genfromtxt(basename+'.positive', dtype=int)
  os.remove(basename+'.positive')
  pulses = pulses.loc[pulses_list]
  if pulses_list.size != pulses.shape[0]:
    raise IndexError('Attention: classified file contains pulses not included!')
  return pulses



def filters(pulses, events, filename, validation=False, header=True):  
  values = pd.DataFrame(dtype=np.float16)
  idx = 0

  events.sort('DM',inplace=True)
  gb = events.groupby('Pulse',sort=False)
  pulses.sort_index(inplace=True)

  values[idx] = (gb.Duration.max() / pulses.Duration).astype(np.float16)
  idx += 1

  def flat_SNR_extremes(sigma):                                            
    dim = np.max((1,sigma.shape[0]/6))
    return np.max((np.median(sigma.iloc[:dim]),np.median(sigma.iloc[-dim:]))) / sigma.max()
  values[idx] = (gb.apply(lambda x: flat_SNR_extremes(x.Sigma))).astype(np.float16)
  idx += 1

  def fit_simm(x,y):
    lim = y.argmax()
    xl = x.loc[:lim]
    if xl.shape[0] < 2: return 10000
    pl = np.polyfit(xl, y.loc[:lim], 1)[0]
    xr = x.loc[lim:]
    if xr.shape[0] < 2: return 10000
    pr = np.polyfit(xr, y.loc[lim:], 1)[0]
    return pl*pr
  values[idx] = (gb.apply(lambda x: fit_simm(x.DM, x.Sigma))).astype(np.float16)
  idx += 1

  #Remove flat duration pulses. Minimum ratio to have weakest pulses with SNR = 8 (from Eq.6.21 of Pulsar Handbook)
  values[idx] = gb.Duration.max() / pulses.Duration - (pulses.Sigma / gb.Sigma.min())**2
  idx += 1

  def extreme_min(ev):
    ev_len = ev.shape[0] / 5
    return np.max((ev[:ev_len].min(), ev[-ev_len:].min()))
  values[idx] = (gb.apply(lambda x: extreme_min(x.Sigma))).astype(np.float16)
  idx += 1

  values[idx] = (gb.Sigma.min() / pulses.Sigma).astype(np.float16)
  idx += 1

  values[idx] = (pulses.dTime / pulses.dDM).astype(np.float16)
  idx += 1

  def std_time(x):
    return np.std(x - x.shift(1))
  values[idx] = (gb.apply(lambda x: std_time(x.Time))).astype(np.float16)
  idx += 1

  values[idx] = (pulses.Sigma - gb.Sigma.min()).astype(np.float16)
  idx += 1

  if validation: values[idx] = (pulses.Pulsar != 'RFI').astype(np.int)
  else: values[idx] = '?%' + np.array(values.index.astype(str))
  
  if header:
    features_list = ''
    for i in range(idx): features_list += '@attribute Feature{} numeric\n'.format(i)
    header = """@relation Training_v3
{}
@attribute class {{0,1}}
@data
    """.format(features_list[:-1])
    with open(filename, 'w') as f:
      f.write(header)
  
  values.to_csv(filename, sep=',', float_format='%10.5f', header=False, index=False, mode='a')

  return





def multimoment(pulses,idL,inc=12):
  pulses.sort(['SAP','BEAM'],inplace=True)
  last_beam = -1
  last_sap = -1
  freq = np.linspace(F_MIN,F_MAX,2592)
  v = 0
  multimoment = np.zeros(pulses.shape[0],dtype=np.float32)
  
  for i,(idx,puls) in enumerate(pulses.iterrows()):
    if (puls.BEAM != last_beam) | (puls.SAP != last_sap):
      beam = puls.BEAM.astype(int)
      sap = puls.SAP.astype(int)
      
      #Open the fits file
      if beam==inc: stokes = 'incoherentstokes'
      else: stokes = 'stokes'
      filename = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
      try: fits = pyfits.open(filename,memmap=True)
      except NameError: 
        logging.warning("RFIexcision - Additional modules missing")
        return
      except IOError: continue
      
      last_beam = beam
      last_sap = sap
  
      header = Utilities.read_header(filename)
      MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
      try: v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
      except NameError: 
        logging.warning("RFIexcision - Additional modules missing")
        return
      
      if puls.DM>141.71: sample = puls.Sample * 4
      elif puls.DM>40.47: sample = puls.Sample * 2
      else: sample = puls.Sample      
      
      sample += np.round(sample*v).astype(int)
      duration = np.int(np.round(puls.Duration/RES))
    
      #Load the spectrum
      spectrum = Utilities.read_fits(fits,puls.DM.copy(),sample.copy(),duration,0)
    
      #De-dispersion
      time = (4149 * puls.DM * (np.power(freq,-2) - F_MAX**-2) / RES).round().astype(np.int)
      spectrum = np.sum([spectrum[time+x,np.arange(2592)] for x in range(duration)],axis=0)
      
      I1 = spectrum.size * np.sum(spectrum**2)
      I2 = np.sum(spectrum)**2
      
      multimoment[i] = (I1 - I2) / I2

  pulses['multimoment'] = multimoment # <= MULTIMOMENT

  return



#CONTROLLARE
beams = {
 13: [14, 15, 16, 17, 18, 19],
 14: [20, 21, 15, 13, 19, 31],
 15: [21, 22, 23, 16, 13, 14],
 16: [15, 23, 24, 25, 17, 13],
 17: [13, 16, 25, 26, 27, 18],
 18: [19, 13, 17, 27, 28, 29],
 19: [31, 14, 13, 18, 29, 30],
 20: [32, 33, 21, 14, 31, 49],
 21: [33, 34, 22, 15, 14, 20],
 22: [34, 35, 36, 23, 15, 21],
 23: [22, 36, 37, 24, 16, 15],
 24: [23, 37, 38, 29, 25, 16],
 25: [16, 24, 39, 40, 26, 17],
 26: [17, 25, 40, 41, 42, 27],
 27: [18, 17, 26, 42, 43, 28],
 28: [29, 18, 27, 43, 44, 45],
 29: [30, 19, 18, 28, 45, 46],
 30: [48, 31, 19, 29, 46, 47],
 31: [49, 20, 14, 19, 30, 48],
 32: [50, 51, 33, 20, 49, 73],
 33: [51, 52, 34, 21, 20, 32],
 34: [52, 53, 35, 22, 21, 33],
 35: [53, 54, 55, 36, 22, 34],
 36: [35, 55, 56, 37, 23, 22],
 37: [36, 56, 57, 38, 24, 23],
 38: [37, 57, 58, 59, 39, 24],
 39: [24, 38, 59, 60, 40, 25],
 40: [25, 39, 60, 61, 41, 26],
 41: [26, 40, 61, 62, 63, 42],
 42: [27, 26, 41, 63, 64, 43],
 43: [28, 27, 42, 64, 65, 44],
 44: [45, 28, 43, 65, 66, 67],
 45: [46, 29, 28, 44, 67, 68],
 46: [47, 30, 29, 45, 68, 69],
 47: [71, 48, 30, 46, 69, 70],
 48: [72, 49, 31, 30, 47, 71],
 49: [73, 32, 20, 31, 48, 72],
 50: [51, 32, 73],
 51: [52, 33, 32, 50],
 52: [53, 34, 33, 51],
 53: [54, 35, 34, 52],
 54: [55, 35, 53],
 55: [56, 36, 35, 54],
 56: [57, 37, 36, 55],
 57: [58, 38, 37, 56],
 58: [59, 38, 57],
 59: [60, 39, 38, 58],
 60: [61, 40, 39, 59],
 61: [62, 41, 40, 60],
 62: [63, 41, 61],
 63: [64, 42, 41, 62],
 64: [65, 43, 42, 63],
 65: [66, 44, 43, 64],
 66: [67, 44, 65],
 67: [68, 45, 44, 66],
 68: [69, 46, 45, 67],
 69: [70, 47, 46, 68],
 70: [71, 47, 69],
 71: [72, 48, 47, 70],
 72: [73, 49, 48, 71],
 73: [50, 32, 49, 72]
}



def time_span(pulses):
  RFI = pd.DataFrame()
  
  for sap in pulses.SAP.unique():
    puls = pulses[pulses.SAP==sap]
      
    try:
      puls_time = puls.Time.round(-1).astype(int)
      puls_time = puls.groupby(puls_time,sort=False)['N_events'].size()  
      mean = puls_time.sum()/360.
      k = stats.poisson.ppf(0.99,mean)
      puls_time = puls_time.index[puls_time>k]
      puls_time = puls.loc[puls.Time.round(-1).astype(int).isin(puls_time),['DM','Time']] 
    except KeyError,AssertionError: puls_time = pd.DataFrame()
    RFI = RFI.append(puls_time)
    
    try:
      puls_time = (puls.Time+5).round(-1).astype(int)
      puls_time = puls.groupby(puls_time,sort=False)['N_events'].size()  
      mean = puls_time.sum()/360.
      k = stats.poisson.ppf(0.99,mean)
      puls_time = puls_time.index[puls_time>k]
      puls_time = puls.loc[puls.Time.round(-1).astype(int).isin(puls_time),['DM','Time']] 
    except KeyError,AssertionError: puls_time = pd.DataFrame()
    RFI = RFI.append(puls_time)
  
  if RFI.empty: return RFI.index
  RFI = RFI.drop_duplicates()
  RFI.sort('Time',inplace=True)
  no_rfi = np.zeros(RFI.shape[0],dtype=np.int8)
  C_Funct.time_span(RFI.DM.astype(np.float32).values,RFI.Time.astype(np.float32).values,no_rfi)
  
  RFI.sort_index(inplace=True)
  return RFI.index[no_rfi==0]


def beam_comparison(pulses, inc=12, database='SinglePulses.hdf5'):
  condition_list_A = [('SAP == sap'), ('BEAM != beam'), ('Time > tmin'), ('Time < tmax'), ('DM > DMmin'), ('DM < DMmax')]
  condition_list_B = '(Sigma >= @SNRmin) & (BEAM != @beams[@beam])'
  
  def comparison(puls, inc):
    #if puls.BEAM <= inc: return 0
    
    sap = int(puls.SAP)
    beam = int(puls.BEAM)
    tmin = float(puls.Time - 2. * puls.Duration)
    tmax = float(puls.Time + 2. * puls.Duration)
    DMmin = float(puls.DM - 0.2)
    DMmax = float(puls.DM + 0.2)
    SNRmin = puls.Sigma / 2.
    
    if pd.read_hdf(database, 'events', where=condition_list_A).query(condition_list_B).groupby('BEAM').count().shape[0] > 3: return 1
    else: return 0

  values = pulses.apply(lambda x: comparison(x, inc), axis=1)
  pulses = pulses.loc[values.index[values == 0]]
  return pulses


