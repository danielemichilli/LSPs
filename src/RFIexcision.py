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
  return pulses



def filters(pulses, events, filename):  
  values = pd.DataFrame(dtype=np.float16)
  idx = 0

  events.sort('DM',inplace=True)
  gb = events.groupby('Pulse',sort=False)
  pulses.sort_index(inplace=True)

  values[idx] = (pulses.Sigma - gb.Sigma.min()).astype(np.float16)
  idx += 1
  
  values[idx] = (gb.Sigma.min() / pulses.Sigma).astype(np.float16)
  idx += 1

  values[idx] = (pulses.dTime / pulses.dDM).astype(np.float16)
  idx += 1

  values[idx] = (np.fabs(pulses.DM-pulses.DM_c)/pulses.dDM).astype(np.float16)
  idx += 1

  #steps = pd.Series()
  #steps = steps.reindex_like(pulses).fillna(0.01)
  #steps[pulses.DM>40.48] = 0.05
  #steps[pulses.DM>141.68] = 0.1
  #values[idx] = (pulses.dDM / (steps * (pulses.N_events - 1))).astype(np.float16)
  #idx += 1

  values[idx] = (gb.Duration.max() / pulses.Duration).astype(np.float16)
  idx += 1

  DM_extremes = pd.DataFrame()
  DM_extremes['Sigma_min'] = gb.Sigma.first()
  DM_extremes['Sigma_max'] = gb.Sigma.last()
  DM_extremes_max = DM_extremes.max(axis=1)
  values[idx] = (DM_extremes_max / pulses.Sigma).astype(np.float16)
  idx += 1

  #def monotonic(y):
    #sigma = np.convolve(y, np.ones(y.shape[0]/5), mode='same')/y.shape[0]*5
    #sigma_max = sigma.argmax()
    #l = sigma[:sigma_max].size*2/3
    #r = sigma[sigma_max:].size*2/3
    #if (l == 0) | (r == 0): return 0
    #sigma = sigma[l:-r]
    #sigma_max = sigma.argmax()
    #sigma = np.diff(sigma)
    #sigma[sigma_max:] *= -1
    #if sigma.size == 1: return sigma[0]
    #return np.partition(sigma,1)[1]
  #values[idx] = (gb.apply(lambda x: monotonic(x.Sigma))).astype(np.float16)
  #idx += 1

  #def sigma_jumps(ev_sigma):
    #sigma = np.convolve(ev_sigma, np.ones(5), mode='same')/5.
    #sigma_max = sigma.argmax()
    #sigma = np.diff(sigma)
    #sigma[sigma_max:] *= -1
    #return sigma[sigma<0].size/float(sigma.size)
  #values[idx] = (gb.apply(lambda x: sigma_jumps(x.Sigma))).astype(np.float16)
  #idx += 1

  #def fit1_brightest(ev):
    #sigma = np.convolve(ev.Sigma, np.ones(3), mode='valid')/3
    #dm = ev.DM.iloc[3/2:-int(3-1.5)/2]
    #sigma = pd.Series(sigma,index=dm.index)
    #DM_c = dm.loc[sigma.argmax()]
    #l = sigma[dm<=DM_c]
    #r = sigma[dm>=DM_c]
    #lim_l = l.min() + np.min((2.,(l.max()-l.min())/4))
    #lim_r = r.min() + np.min((2.,(r.max()-r.min())/4))
    #l = l[l>lim_l]
    #r = r[r>lim_r]
    #y = pd.concat((l,r))
    #x = dm.loc[y.index]
    #if (y.unique().size < 2) & (x.unique().size < 2): return 0
    #p = np.polyfit(x, y, 1)
    #return np.sum((np.polyval(p, x) - y) ** 2) / (x.size-1)
  #values[idx] = (gb.apply(lambda x: fit1_brightest(x))).astype(np.float16)
  #idx += 1

  def extreme_min(ev):
    ev_len = ev.shape[0] / 4
    SNR_min_a = ev[:ev_len].min()
    SNR_min_b = ev[-ev_len:].min()
    return np.max((SNR_min_a, SNR_min_b))
  values[idx] = (gb.apply(lambda x: extreme_min(x.Sigma))).astype(np.float16)
  idx += 1

  def mean_SNR(x):
    return np.mean(np.fabs( x - x.shift(-1) ) / x)
  values[idx] = (gb.apply(lambda x: mean_SNR(x.Sigma))).astype(np.float16)
  idx += 1

  def mean_time(x):
    return np.mean(np.abs(x - x.shift(1)))
  values[idx] = (gb.apply(lambda x: mean_time(x.Time))).astype(np.float16)
  idx += 1

  #def std_time(x):
    #return np.std(x - x.shift(1))
  #values[idx] = (gb.apply(lambda x: std_time(x.Time))).astype(np.float16)
  #idx += 1

  #def sigma_std_largest(ev):
    #sigma = ev.Sigma.nlargest(ev.Sigma.size*2/3)
    #return np.std(sigma)
  #values[idx] = (gb.apply(lambda x: sigma_std_largest(x))).astype(np.float16)
  #idx += 1

  #def fit0(x,y):
    #p = np.polyfit(x, y, 0)
    #return np.sum((np.polyval(p, x) - y) ** 2) / x.size
  #values[idx] = (gb.apply(lambda x: fit0(x.DM, x.Sigma))).astype(np.float16)
  #idx += 1

  #def fit1(x,y):
    #p = np.polyfit(x, y, 1)
    #return np.sum((np.polyval(p, x) - y) ** 2) / x.size
  #values[idx] = (gb.apply(lambda x: fit1(x.DM, x.Sigma))).astype(np.float16)
  #idx += 1

  def SNR_simmetric(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    r = ev[ev.DM>=DM_c]
    return np.max((l.Sigma.min(),r.Sigma.min())) / ev.Sigma.max()
  values[idx] = (gb.apply(lambda x: SNR_simmetric(x))).astype(np.float16)
  idx += 1

  def bright_events_abs(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    r = ev[ev.DM>=DM_c]
    lim_l = l.Sigma.min() + np.min((2.,(l.Sigma.max()-l.Sigma.min())/4))
    lim_r = r.Sigma.min() + np.min((2.,(r.Sigma.max()-r.Sigma.min())/4))
    l = l[l.Sigma>lim_l]
    r = r[r.Sigma>lim_r]
    ev = pd.concat((l,r))
    ev.drop_duplicates(inplace=True)
    return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()
  values[idx] = (gb.apply(lambda x: bright_events_abs(x))).astype(np.float16)
  idx += 1

  #def bright_events_rel(ev):
    #DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    #l = ev[ev.DM<=DM_c]
    #if l.shape[0] > 0:
      #l_lim = np.cumsum(l.Sigma-l.Sigma.iloc[0])
      #l = l[l_lim >= ev.Sigma.max()/8.]
    #r = ev[ev.DM>DM_c]
    #if r.shape[0] > 0:
      #r.sort('DM',inplace=True,ascending=False)
      #r_lim = np.cumsum(r.Sigma-r.Sigma.iloc[0])
      #r = r[r_lim >= ev.Sigma.max()/8.]
    #ev = pd.concat((l,r))
    #ev.drop_duplicates(inplace=True)
    #if ev.shape[0] == 0: return 1
    #return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()
  #values[idx] = (gb.apply(lambda x: bright_events_rel(x))).astype(np.float16)
  #idx += 1

  #def pulse_simmetric(ev):
    #DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    #x = ev.DM[ev.DM<=DM_c]
    #if x.size < 2: return 0
    #y = ev.Sigma[ev.DM<=DM_c]
    #ml = np.polyfit(x, y, 1)[0]
    #if ml < 1e-7: return 0
    #x = ev.DM[ev.DM>=DM_c]
    #if x.size < 2: return 0
    #y = ev.Sigma[ev.DM>=DM_c]
    #mr = np.polyfit(x, y, 1)[0]
    #if mr < 1e-7: return 0
    #return np.min((-ml/mr,-mr/ml))
  #values[idx] = (gb.apply(lambda x: pulse_simmetric(x))).astype(np.float16)
  #idx += 1

  def flat_SNR_extremes(ev):                                            
    return np.max((ev.Sigma.iloc[1],ev.Sigma.iloc[-2]))/ev.Sigma.max()
  values[idx] = (gb.apply(lambda x: flat_SNR_extremes(x))).astype(np.float16)
  idx += 1

  def number_events(ev):
    dim = ev.shape[0]/5
    sigma = np.convolve(ev.Sigma, np.ones(dim), mode='valid')/dim
    dm = ev.DM.iloc[dim/2:np.min((-int(dim-1.5)/2,-1))]
    sigma_argmax = sigma.argmax()
    sigma_max = sigma.max()
    try: lim_max = np.max((sigma[:sigma_argmax].min(),sigma[sigma_argmax:].min()))
    except ValueError: return 0
    lim_max = lim_max+(sigma_max-lim_max)/5.
    l = np.where(sigma[:sigma_argmax]<=lim_max)[0][-1]+1
    r = (np.where(sigma[sigma_argmax:]<=lim_max)[0]+sigma_argmax)[0]-1
    duration = np.convolve(ev.Duration, np.ones(dim), mode='valid')/dim
    duration = duration[sigma_argmax]
    try: dDM = dm.iloc[sigma_argmax] - dm.iloc[l]
    except IndexError: return 0
    y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
    diff_l = lim_max/sigma_max/y
    dDM = dm.iloc[r] - dm.iloc[sigma_argmax]
    y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
    diff_r = lim_max/sigma_max/y
    if np.isnan(diff_l) & np.isnan(diff_r): return 0
    return np.nanmax((diff_l,diff_r))
  values[idx] = (gb.apply(lambda x: number_events(x))).astype(np.float16)
  idx += 1

  def crosses(sig):
    diff = sig - (sig.max() + sig.min()) / 2.
    count = np.count_nonzero(np.diff(np.sign(diff)))
    return count
  values[idx] = (gb.apply(lambda x: crosses(x.Sigma))).astype(np.float16)
  idx += 1

  values[idx] = (values[idx-1] % 2).astype(np.float16)
  idx += 1

  #DM_extremes_min = DM_extremes.min(axis=1)
  #values[idx] = (gb.Sigma.min() / DM_extremes_min).astype(np.float16)
  #idx += 1

  values[idx] = '?%' + np.array(values.index.astype(str))


  features_list = ''
  for i in range(idx): features_list += '@attribute Feature{} numeric\n'.format(i)
  header = """@relation Training_v2
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


#def time_span(puls):
  #def min_puls(x):   
    #if x.size <= 1: return 0
    #else: return np.min(np.abs(x-np.mean(x)))

  #lim = 1-0.01/3600*puls.shape[0]
  
  #try:
    #puls_time = puls.Time.astype(int)
    #puls_time = puls.groupby(['SAP',puls_time]).agg({'N_events':np.size,'DM':min_puls})  
    #mean = puls_time.N_events.sum()/3600.
    #k = stats.poisson.ppf(lim,mean)   
    #puls_time = puls_time.index[(puls_time.N_events>k)&(puls_time.DM>1)].get_level_values('Time')
    #a = puls.Pulse[puls.Time.astype(int).isin(puls_time)] 
  #except KeyError,AssertionError: a = pd.DataFrame()
  
  #try:
    #puls_time = (puls.Time+0.5).astype(int)
    #puls_time = puls.groupby(['SAP',puls_time]).agg({'N_events':np.size,'DM':min_puls})
    #mean = puls_time.N_events.sum()/3600.
    #k = stats.poisson.ppf(lim,mean)
    #puls_time = puls_time.index[(puls_time.N_events>k)&(puls_time.DM>1)].get_level_values('Time')
    #b = puls.Pulse[(puls.Time+0.5).astype(int).isin(puls_time)] 
  #except KeyError,AssertionError: b = pd.DataFrame()
  
  #puls_time = pd.concat((a,b)).index.unique()
  #puls_time.sort()
    
  #return puls_time


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


def puls_beams_select(x,puls):
  if x.BEAM<13: return 0
  select = puls.BEAM[(puls.SAP==x.SAP)&(puls.BEAM!=x.BEAM)&(np.abs(puls.Time-x.Time)<=0.1)&(np.abs(puls.DM-x.DM)<=1.)]
  close = select[select.isin(beams[x.BEAM])].size
  away = select.size - close
  if x.Sigma <= 13: return np.max((close>3,away>1))
  elif x.Sigma <= 19:  #for sidelobe sensitivity 0.3 of main lobe 
    return away>1
  else: return 0


def Compare_Beams(puls):
  puls.Pulse[:] = 0
    
  sap0 = puls[puls.SAP==0].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap0['Time_low'] = sap0.Time_c-sap0.dTime
  sap0.sort('Time_low',inplace=True)
  
  sap1 = puls[puls.SAP==1].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap1['Time_low'] = sap1.Time_c-sap1.dTime
  sap1.sort('Time_low',inplace=True)
  
  sap2 = puls[puls.SAP==2].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap2['Time_low'] = sap2.Time_c-sap2.dTime
  sap2.sort('Time_low',inplace=True)
  
  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,np.int8(1))
  
  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  C_Funct.Compare(sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  #puls.Pulse.loc[puls.SAP==0] = sap0.Pulse
  
  #puls.Pulse.loc[puls.SAP==1] = sap1.Pulse
  
  #puls.Pulse.loc[puls.SAP==2] = sap2.Pulse
  
  idx = pd.concat((sap0[sap0.Pulse>0],sap1[sap1.Pulse>0],sap2[sap2.Pulse>0]))
  idx.sort_index(inplace=True)
  
  return idx.index


