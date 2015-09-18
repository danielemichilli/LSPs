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
from scipy.optimize import curve_fit
from scipy import special

import C_Funct
from Parameters import *

import time

def Pulse_Thresh(pulses,events):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  events.DM = events.DM.astype(np.float64)
  events.Sigma = events.Sigma.astype(np.float64)
  
  events.sort_index(inplace=True)
  gb = events.groupby('Pulse')
  pulses.sort_index(inplace=True)
  
  #Weak pulses
  pulses.Pulse[pulses.Sigma - gb.Sigma.min() <= 1.5] += RFI_percent
  
  #Pulses peaked at low DM
  pulses.Pulse[pulses.DM <= 3.2] += RFI_percent

  #Time-aligned: dTime / pulses.dDM
  pulses.Pulse[pulses.dTime / pulses.dDM > FILTERS['aligned']] += 1

  #Peak central: | DM - DM_c | / dDM
  pulses.Pulse[np.fabs(pulses.DM-pulses.DM_c)/pulses.dDM > FILTERS['peak_central']] += 1

  #Missing events: 2 * dDM / (bin * (N - 1))
  steps = pd.Series()
  steps = steps.reindex_like(pulses).fillna(0.01)
  steps[pulses.DM>40.48] = 0.05
  steps[pulses.DM>141.68] = 0.1
  pulses.Pulse[pulses.dDM / (steps * (pulses.N_events - 1)) > FILTERS['holes']] += 1
  steps = 0

  #Flat Duration: Duration_min / Duration_max
  pulses.Pulse[ gb.Duration.min() / gb.Duration.max() > FILTERS['flat_duration']] += 1 

  #Flat SNR: SNR_min / SNR
  pulses.Pulse[gb.Sigma.min() / pulses.Sigma > FILTERS['flat_SNR']] += 1   #Da testare meglio: alcune volte elimina i pulse
  
  #Strong extreme events
  DM_extremes = pd.DataFrame()
  DM_extremes['Sigma_min'] = gb.Sigma.first()
  DM_extremes['Sigma_max'] = gb.Sigma.last()
  DM_extremes_max = DM_extremes.max(axis=1)
  pulses.Pulse[ DM_extremes_max / pulses.Sigma > FILTERS['DM_extremes']] += 1
  DM_extremes_max = 0
  
  #Minimum different from extremes
  DM_extremes_min = DM_extremes.min(axis=1)
  pulses.Pulse[ gb.Sigma.min() / DM_extremes_min <  FILTERS['sigma_min']] += 1  #forse si puo' implementare con diverse soglie sulla sigma
  DM_extremes_min = 0
  DM_extremes = 0
  
  def pulses_apply(ev):
    s1 = ev.Sigma - ev.Sigma.shift(-1)
    s2 = ev.Sigma - ev.Sigma.shift(1)
    s1.fillna(0,inplace=True)
    s2.fillna(0,inplace=True)
    s = pd.concat((s1[s1<s2],s2[s2<=s1])).sort_index()
    ev = ev[s>-5]
    if ev.shape[0] < 5: return 1
    ev.sort('DM',inplace=True)  #controllare
    return np.sum((
      np.mean(np.fabs( ev.Sigma - ev.Sigma.shift(-1) ) / ev.Sigma) < FILTERS['sigma_scatter'],
      (np.mean(np.abs(ev.Time - ev.Time.shift(1))) > FILTERS['cum_scatter']) |
      (np.std(ev.Time - ev.Time.shift(1)) > FILTERS['std_scatter']),
      sigma_std_largest(ev) | fit0(ev.DM,ev.Sigma) | fit1(ev.DM,ev.Sigma),
      SNR_simmetric(ev) / ev.Sigma.max() > FILTERS['flat_SNR_simmetric'],
      bright_events_abs(ev) > FILTERS['bright_extremes_abs'],
      bright_events_rel(ev) > FILTERS['bright_extremes_rel'],
      pulse_simmetric(ev) < FILTERS['pulse_simmetric'],
      flat_SNR_extremes < FILTERS['flat_SNR_extremes'],
      number_events(ev) < FILTERS['number_events'],
      monotonic(ev.Sigma) < FILTERS['monotonic'],
      sigma_jumps(ev.Sigma) > FILTERS['sigma_jumps'],
      fit1_brightest(ev) < FILTERS['fit1_brightest']))
 
  def sigma_std_largest(ev):
    sigma = ev.Sigma.nlargest(ev.Sigma.size*2/3)
    if sigma.size<20: return 0
    if sigma.max()<8: return np.std(sigma) < FILTERS['sigma_std_largest_weak']
    else: return np.std(sigma) < FILTERS['sigma_std_largest']
  

  def fit0(x,y):
    if x.size<20: return 0
    p = np.polyfit(x, y, 0)
    if y.max()<8: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit0_weak']
    else: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit0']
  
  def fit1(x,y):
    if x.size<20: return 0
    p = np.polyfit(x, y, 1)
    if y.max()<8: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit1_weak']
    else: return np.sum((np.polyval(p, x) - y) ** 2) / x.size < FILTERS['flat_fit1']
  
  
  def pulse_simmetric(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    x = ev.DM[ev.DM<=DM_c]
    y = ev.Sigma[ev.DM<=DM_c]
    ml = np.polyfit(x, y, 1)[0]
    x = ev.DM[ev.DM>=DM_c]
    y = ev.Sigma[ev.DM>=DM_c]
    mr = np.polyfit(x, y, 1)[0]
    return np.min((-ml/mr,-mr/ml))
  
  
  def number_events(ev):
    dim = ev.shape[0]/5
    if dim < 3: return 10
    sigma = np.convolve(ev.Sigma, np.ones(dim), mode='valid')/dim
    dm = ev.DM.iloc[dim/2:-int(dim-1.5)/2]
    sigma_argmax = sigma.argmax()
    sigma_max = sigma.max()
    try: lim_max = np.max((sigma[:sigma_argmax].min(),sigma[sigma_argmax:].min()))
    except ValueError: return 0
    lim_max = lim_max+(sigma_max-lim_max)/5.
    l = np.where(sigma[:sigma_argmax]<=lim_max)[0][-1]+1
    r = (np.where(sigma[sigma_argmax:]<=lim_max)[0]+sigma_argmax)[0]-1
    if (sigma_argmax - l < 5) & (r - sigma_argmax < 5): return 10
    duration = np.convolve(ev.Duration, np.ones(dim), mode='valid')/dim
    duration = duration[sigma_argmax]
    dDM = dm.iloc[sigma_argmax] - dm.iloc[l]
    y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
    diff_l = lim_max/sigma_max/y
    dDM = dm.iloc[r] - dm.iloc[sigma_argmax]
    y = np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/duration/0.13525**3)*special.erf(0.00000691*dDM*31.64/duration/0.13525**3)
    diff_r = lim_max/sigma_max/y
    return np.nanmax((diff_l,diff_r))

  def monotonic(y):
    sigma = np.convolve(y, np.ones(y.shape[0]/5), mode='same')/y.shape[0]*5
    sigma_max = sigma.argmax()
    l = sigma[:sigma_max].size*2/3
    r = sigma[sigma_max:].size*2/3
    sigma = sigma[l:-r]
    if sigma.size < 10: return 1
    sigma_max = sigma.argmax()
    sigma = np.diff(sigma)
    sigma[sigma_max:] *= -1
    return np.partition(sigma,1)[1]


  def sigma_jumps(ev_sigma):
    sigma = np.convolve(ev_sigma, np.ones(5), mode='same')/5.
    sigma_max = sigma.argmax()
    sigma = np.diff(sigma)
    sigma[sigma_max:] *= -1
    return sigma[sigma<0].size/float(sigma.size)


  def SNR_simmetric(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    r = ev[ev.DM>=DM_c]
    return np.max((l.Sigma.min(),r.Sigma.min()))
  
  
  def window_sum_SNR(Sigma):
    val = np.convolve(Sigma,np.ones(3),mode='valid')
    return val.min()/val.max()
   

  def fit1_brightest(ev):
    sigma = np.convolve(ev.Sigma, np.ones(3), mode='valid')/3
    dm = ev.DM.iloc[3/2:-int(3-1.5)/2]
    sigma = pd.Series(sigma,index=dm.index)
    DM_c = dm.loc[sigma.argmax()]
    l = sigma[dm<=DM_c]
    if l.size<=4: return 10
    r = sigma[dm>=DM_c]
    if r.size<4: return 10
    lim_l = l.min() + np.min((2.,(l.max()-l.min())/4))
    lim_r = r.min() + np.min((2.,(r.max()-r.min())/4))
    l = l[l>lim_l]
    r = r[r>lim_r]
    y = pd.concat((l,r))
    if y.size<=5: return 10
    x = dm.loc[y.index]
    #if y.max()<8: return 10
    p = np.polyfit(x, y, 1)
    return np.sum((np.polyval(p, x) - y) ** 2) / (x.size-1)
    
  
  #rimuove gli eventi piu' deboli a destra e sinistra.
  def bright_events_abs(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    if l.shape[0]<=4: return 0
    r = ev[ev.DM>=DM_c]
    if r.shape[0]<4: return 0
    lim_l = l.Sigma.min() + np.min((2.,(l.Sigma.max()-l.Sigma.min())/4))
    lim_r = r.Sigma.min() + np.min((2.,(r.Sigma.max()-r.Sigma.min())/4))
    l = l[l.Sigma>lim_l]
    r = r[r.Sigma>lim_r]
    ev = pd.concat((l,r))
    if ev.shape[0]<=5: return 0
    try: return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()
    except ValueError: return 1

  def bright_events_rel(ev):
    if ev.Sigma.max()<8: return 0
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    if l.shape[0]<=4: return 0
    r = ev[ev.DM>DM_c]
    if r.shape[0]<4: return 0
    r.sort('DM',inplace=True,ascending=False)
    l_lim = np.cumsum(l.Sigma-l.Sigma.iloc[0])
    r_lim = np.cumsum(r.Sigma-r.Sigma.iloc[0])
    l = l[l_lim >= ev.Sigma.max()/8.]
    r = r[r_lim >= ev.Sigma.max()/8.]
    ev = pd.concat((l,r))
    if ev.shape[0]<5: return 0
    else: return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()


  def flat_SNR_extremes(ev):                                            
    if ev.shape[0] < 30: return 0
    else: return np.max((ev.Sigma.iloc[1],ev.Sigma.iloc[-2]))/ev.Sigma.max()
  
  
  #pulses = pulses[pulses.Pulse<=RFI_percent]
  #events = events[events.Pulse.isin(pulses.index)]
  #gb = events.groupby('Pulse')
  pulses.Pulse += gb.apply(lambda x: pulses_apply(x)).astype(np.int8)
  
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


def time_span(puls):
  #TESTARE
  def min_puls(x):   
    if x.size <= 1: return 0
    else: return np.min(np.abs(x-np.mean(x)))
  
  try:
    puls_time = puls.Time.astype(int)
    puls_time = puls.groupby(['SAP',puls_time]).agg({'N_events':np.size,'DM':min_puls})
    puls_time = puls_time.index[(puls_time.N_events>=5)&(puls_time.DM>1)].get_level_values('Time')  #va bene 5??
    a = puls.Pulse[puls.Time.astype(int).isin(puls_time)] 
  except AssertionError: a = pd.DataFrame()
  
  try:
    puls_time = (puls.Time+0.5).astype(int)
    puls_time = puls.groupby(['SAP',puls_time]).agg({'N_events':np.size,'DM':min_puls})
    puls_time = puls_time.index[(puls_time.N_events>=5)&(puls_time.DM>1)].get_level_values('Time')
    b = puls.Pulse[(puls.Time+0.5).astype(int).isin(puls_time)] 
  except AssertionError: b = pd.DataFrame()
  
  puls_time = pd.concat((a,b)).index.unique()
  puls.Pulse.loc[puls_time] += 1
  
  return




def Compare_Beams(puls):
  
  #STUDIARE VALORI E SE METTERE puls.Pulse==0
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=36000)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=3600)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=360)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  
  
  time_span(puls)
  
  

  def puls_beams_select(x):
    select = puls.BEAM[(puls.SAP==x.SAP)&(puls.BEAM!=x.BEAM)&(np.abs(puls.Time-x.Time)<=0.1)&(np.abs(puls.DM-x.DM)<=1.)]
    close = select[select.isin(beams[x.BEAM])].size
    away = select.size - close
    if x.Sigma <= 13: return np.sum((close>3,away>1))
    elif x.Sigma <= 20: 
      if len(beams[x.BEAM])==6: return np.sum((close<3,away>1))
      else: return np.sum(away>1)
    else: 
      if len(beams[x.BEAM])==6: return np.sum(close<3)
      else: return 0  
  
  puls.Pulse += puls.apply(lambda x: puls_beams_select(x),axis=1).astype(np.int8)

  
  sap0 = puls[puls.SAP==0].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap0['Time_low'] = sap0.Time_c-sap0.dTime
  sap0.sort('Time_low',inplace=True)
  
  sap1 = puls[puls.SAP==1].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap1['Time_low'] = sap1.Time_c-sap1.dTime
  sap1.sort('Time_low',inplace=True)
  
  sap2 = puls[puls.SAP==2].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap2['Time_low'] = sap2.Time_c-sap2.dTime
  sap2.sort('Time_low',inplace=True)
  
  logging.info('Comparison is starting')

  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,np.int8(1))
  
  logging.info('1/3 completed')
  
  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  logging.info('2/3 completed')
  
  C_Funct.Compare(sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  puls.Pulse.loc[puls.SAP==0] = sap0.Pulse
  
  puls.Pulse.loc[puls.SAP==1] = sap1.Pulse
  
  puls.Pulse.loc[puls.SAP==2] = sap2.Pulse
    
  return



def Multimoment(ds,DM,duration):
  freq = np.linspace(F_MIN,F_MAX,2592)
  time = (4149 * DM * (np.power(freq,-2) - F_MAX**-2) / RES + DS_OFFSET).round().astype(np.int)
  
  spectrum = ds[100:]
  spectrum = np.sum([spectrum[time+x,np.arange(2592)] for x in range(duration)],axis=0)
  
  I1 = spectrum.size * np.sum(spectrum**2)
  I2 = np.sum(spectrum)**2
  
  return (I1 - I2) / I2

