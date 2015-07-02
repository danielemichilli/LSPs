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
  
  events.sort_index(inplace=True)
  gb = events.groupby('Pulse')
  pulses.sort_index(inplace=True)
    
  #Not scattered in time: dTime
  #pulses.Pulse[pulses.dTime > FILTERS['scattered']] += 1  #inutile

  #Time-aligned: dTime / pulses.dDM
  pulses.Pulse[pulses.dTime / pulses.dDM > FILTERS['aligned']] += 1

  #Peak central: | DM - DM_c | / dDM
  pulses.Pulse[np.fabs(pulses.DM-pulses.DM_c)/pulses.dDM > FILTERS['peak_central']] += 1

  #Duration_min central: Duration > median(Duration)
  #pulses.Pulse[pulses.Duration > gb.Duration.mean() * FILTERS['duration_central']] += 1 

  #Missing events: 2 * dDM / (bin * (N - 1))
  steps = pd.Series()
  steps = steps.reindex_like(pulses).fillna(0.01)
  steps[pulses.DM>40.48] = 0.05
  steps[pulses.DM>141.68] = 0.1
  pulses.Pulse[pulses.dDM / (steps * (pulses.N_events - 1)) > FILTERS['holes']] += 1
  steps = 0

  #Variance on time: std( Time )
  #events.Time = events.Time.astype(np.float64)
  #pulses.Pulse[gb.Time.std() > FILTERS['variance']] += 1  #inutile

  #Flat Duration: Duration_min / Duration_max
  pulses.Pulse[ gb.Duration.min() / gb.Duration.max() > FILTERS['flat_duration']] += 1 

  #Large weak pulses: N_events > Sigma *m+q
  pulses.Pulse[pulses.N_events > pulses.Sigma * FILTERS['m'] + FILTERS['q']] += 1
  
  #dDM > Sigma *m+q

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
  
  #Sigma variance
  pulses.Pulse[gb.Sigma.std() < FILTERS['sigma_std']] += 1  
  
  
  
  
  

  def pulses_apply(ev):
    s1 = ev.Sigma - ev.Sigma.shift(-1)
    s2 = ev.Sigma - ev.Sigma.shift(1)
    
    s1.fillna(0,inplace=True)
    s2.fillna(0,inplace=True)
    s = pd.concat((s1[s1<s2],s2[s2<=s1])).sort_index()
    
    #s = np.nanmin((s1,s2),axis=0)
    ev = ev[s>-5]
    
    if ev.shape[0] < 5: return 1
    
    #s = s[s>-5]
    #s = np.abs(s)/ev.Sigma
    
    s = np.fabs( ev.Sigma - ev.Sigma.shift(-1) ) / ev.Sigma
    
    return np.sum((
      #np.mean(s) < FILTERS['sigma_scatter'],
      #np.mean(np.abs(ev.Time - ev.Time.shift(1))) > FILTERS['cum_scatter'], \  #mettere un filtro sullo scattering in tempo tipo a come fatto per la sigma
      np.std(ev.Sigma.nlargest(ev.Sigma.size*2/3)) < FILTERS['sigma_std_largest'], #OK  #Aggiungere altri cosi
      fit0(ev.DM,ev.Sigma) < FILTERS['flat_fit0'],  #OK
      fit1(ev.DM,ev.Sigma) < FILTERS['flat_fit1'],  #OK
      SNR_simmetric(ev) / ev.Sigma.max() > FILTERS['flat_SNR_simmetric'],  #OK
      bright_events_abs(ev) > FILTERS['bright_extremes_abs'],
      bright_events_rel(ev) > FILTERS['bright_extremes_rel'],
      #fit(ev) < FILTERS['chi2'],
      flat_SNR_extremes > FILTERS['flat_SNR_extremes'],
      number_events(ev) > FILTERS['number_events'],
      monotonic(ev.Sigma) < FILTERS['monotonic']))
      
    
  def fit(ev):
    if ev.Sigma.max()>10: return 0
    
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    p0 = [0.,5.,2.,4.]
    
    l = ev[ev.DM<=DM_c]
    if l.shape[0]-len(p0) <= 0 : return 0
    else:
      x1 = l.DM - l.DM.min()
      y1 = l.Sigma - l.Sigma.min()
      try: 
        p1, var = curve_fit(exp, x1, y1, p0=p0, maxfev=10000)
      except RuntimeError: return np.inf
      #x1 = np.linspace(0,x1.max(),10)
      #y1 = exp(x1, *p1)

    r = ev[ev.DM>=DM_c]
    if r.shape[0]-len(p0) <= 0 : return 0
    else:
      x2 = r.DM.max() - r.DM
      y2 = r.Sigma - r.Sigma.min()
      try: 
        p2, var = curve_fit(exp, x2, y2, p0=p0, maxfev=10000)
      except RuntimeError: return np.inf
      #x2 = np.linspace(0,x2.max(),10)
      #y2 = exp(x2, *p2)
      
    chi_l = np.sum((exp(x1, *p1) - y1) ** 2) / x1.shape[0]
    chi_r = np.sum((exp(x2, *p2) - y2) ** 2) / x2.shape[0]
    
    return np.max((chi_l,chi_r))
      

  def exp(x, *p):             
        A,b,C,D = p
        return A*np.exp(b*x)+C*x**2+D*x

  def fit0(x,y):
    p = np.polyfit(x, y, 0)
    chi_squared = np.sum((np.polyval(p, x) - y) ** 2)
    chi = chi_squared / x.size
    return chi

  def fit1(x,y):
    p = np.polyfit(x, y, 1)
    chi_squared = np.sum((np.polyval(p, x) - y) ** 2)
    chi = chi_squared / (x.size-1)
    return chi
  
  
  #non funziona piu'!!!!
  def number_events(ev):  #forse prendere punti interni invece che sopra soglia
    y1 = np.convolve(ev.Sigma, np.ones(10), mode='same')/10.
    #try:
      #l = np.where(y1[:y1.argmax()]<=6)[0][-1]
      #r = (np.where(y1[y1.argmax():]<=6)[0]+y1.argmax())[0]
    #except IndexError: return 0
    
    
    l = y1[:sigma_max].size*2/3  #da provare
    r = y1[sigma_max:].size*2/3    
    
    
    dDM = ev.DM.iloc[r]-ev.DM.iloc[l]
    return y1.max() / 6.*(np.sqrt(np.pi)/2/(0.00000691*dDM*31.64/ev.Duration.min()/0.13525**3)*special.erf(0.00000691*dDM*31.64/ev.Duration.min()/0.13525**3))


  def chi2(ev):
    y1 = np.convolve(ev.Sigma, np.ones(10), mode='same')/10.
    sigma_max = y1.argmax()
    
    l = y1[:sigma_max].size*2/3
    r = y1[sigma_max:].size*2/3
    
    y1 = y1[l:-r]
    if y1.shape[0] < 10: return 0
    #sigma_max = y1.argmax()
    
    dDM = np.abs(ev.DM.iloc[l:-r]-ev.DM.iloc[sigma_max])
    x = 0.00000691*dDM*31.64/ev.Duration.min()/0.13525**3+0.0000001
    y = y1.max()*np.sqrt(np.pi)/2/x*special.erf(x)
    return np.sum((y - y1) ** 2) / y.shape[0]
  
  
  
  
  #def monotonic(x):
    #lunghezza = x.shape[0]/5
    #a = x[:lunghezza].mean().round(1)
    #b = x[lunghezza:2*lunghezza].mean().round(1)
    #c = x[2*lunghezza:-2*lunghezza].mean().round(1)
    #d  = x[-2*lunghezza:-lunghezza].mean().round(1)
    #e  = x[-lunghezza:].mean().round(1)
    #if (c>=d)&(d>=e)&(c>=b)&(b>=a):
        #return False
    #elif (b>=a)&(b>=c)&(c>=d)&(d>=e):
        #return False
    #elif (d>=e)&(d>=c)&(c>=b)&(b>=a):
        #return False
    #else: return True


  def monotonic(sigma):  #testare
    sigma = np.convolve(sigma, np.ones(sigma.shape[0]/5), mode='same')/sigma.shape[0]*5
    sigma_max = sigma.argmax()
    
    l = sigma[:sigma_max].size*2/3
    r = sigma[sigma_max:].size*2/3
    
    sigma = sigma[l:-r]
    if sigma.size < 10: return 1
    
    sigma_max = sigma.argmax()
    sigma = np.diff(sigma)
    sigma[sigma_max:] *= -1
    #return sigma.min()
    return np.partition(sigma,1)[1]



  def SNR_simmetric(ev):
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    r = ev[ev.DM>=DM_c]
    return np.max((l.Sigma.min(),r.Sigma.min()))
  
  
  def window_sum_SNR(Sigma):
    val = np.convolve(Sigma,np.ones(3),mode='valid')
    return val.min()/val.max()
  


  #rimuove gli eventi piu' deboli a destra e sinistra.
  #si possono applicare diverse statistiche, es. fit1 o meglio rapporto tra sigma degli estremi
  def bright_events_abs(ev):  #da testare!
    DM_c = ev.DM.loc[ev.Sigma.idxmax()]
    l = ev[ev.DM<=DM_c]
    if l.shape[0]<=4: return 0
    r = ev[ev.DM>DM_c]
    if r.shape[0]<4: return 0
    l = l[l.Sigma>l.Sigma.min()+1]
    r = r[r.Sigma>r.Sigma.min()+1]
    ev = pd.concat((l,r))
    if ev.shape[0]<5: return 0
    else: return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()


  def bright_events_rel(ev):
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
  pulses.Pulse += gb.apply(lambda x: pulses_apply(x))
  
  
  return

#Meta' degli eventi piu' vicini al centro: chi2 e altre statistiche
#Meta' degli eventi con sigma piu' alta: chi2




def SNR_mean(x):
  ev = x[x.Sigma>x.Sigma.mean()]
  DM_max = ev.DM.max()
  DM_min = ev.DM.min()
  return np.fabs(ev.DM.loc[ev.Sigma.idxmax()]-((DM_max+DM_min)/2) )/(DM_max-DM_min)


def SNR_mean(x):
  ev = x[x.Sigma>=x.Sigma.median()]
  DM_extremes = np.max((ev.Sigma.loc[ev.DM.idxmin()],ev.Sigma.loc[ev.DM.idxmax()]))
  return DM_extremes / ev.Sigma.max()



def sigma_jumps(x):
  x.loc[0] = 5
  snr = x - x.shift(-1)
  snr[np.abs(snr)>0.2*x.max()] = 0
  return np.fabs(snr.sum())


def sigma_jumps(x):  #forse meglio con rapporto
  snr = x / x.shift(-1)
  snr = snr[(snr>1.16)|(snr<1/1.16)]
  prod = np.prod(snr)
  if np.isnan(prod): prod = 0
  return prod
  




#rimuove tanti eventi alla stessa bassa sigma
def low_events(ev):
  DM_c = ev.DM.loc[ev.Sigma.idxmax()]
  l = ev[ev.DM<=DM_c]
  low_l = l[l.Sigma<=l.Sigma.min()+1]
  r = ev[ev.DM>=DM_c]
  low_r = r[r.Sigma<=r.Sigma.min()+1]
  low_l = np.float(low_l.shape[0])/l.shape[0]
  low_r = np.float(low_r.shape[0])/r.shape[0]
  return np.max((low_l,low_r))



def central_events(ev):
  DM_c = ev.DM.loc[ev.Sigma.idxmax()]
  l = ev[ev.DM<=DM_c]
  r = ev[ev.DM>DM_c]
  l = l[l.DM>=(l.DM.min()+l.DM.max())/2]
  r = r[r.DM>=(r.DM.min()+r.DM.max())/2]
  ev = pd.concat((l,r))
  if ev.empty: return 0
  else: return np.max((ev.Sigma[ev.DM.argmin()],ev.Sigma[ev.DM.argmax()]))/ev.Sigma.max()



def simmetric(ev):
  DM_c = ev.DM.loc[ev.Sigma.idxmax()]
  l = ev[ev.DM<DM_c]
  r = ev[ev.DM>DM_c]
  return np.max((np.float(l.shape[0])/r.shape[0],np.float(r.shape[0])/l.shape[0]))





##prove
##Expected shape
#dDM = np.abs(ev.DM-ev.DM.loc[ev.Sigma.argmax()])
#x = 0.00000691*dDM*31.64/ev.Duration.min()/0.13525**3+0.00001
#y = ev.Sigma.max()*np.sqrt(np.pi)/2/x*special.erf(x)

##fit
#y1 = np.convolve(ev.Sigma, np.ones(10), mode='same')/10.
#dDM = np.abs(ev.DM-ev.DM.iloc[y1.argmax()])
#x = 0.00000691*dDM*31.64/ev.Duration.iloc[y1.argmax()]/0.13525**3+0.00001
#y = y1.max()*np.sqrt(np.pi)/2/x*special.erf(x)

##number
#def number(ev):
  #y1 = np.convolve(ev.Sigma, np.ones(10), mode='same')/10.
  #l = np.where(y1[:y1.argmax()]<=6)[0][-1]
  #r = (np.where(y1[y1.argmax():]<=6)[0]+y1.argmax())[0]
  #dDM = ev.DM.iloc[r]-ev.DM.iloc[l]
  #x = 0.00000691*dDM*31.64/ev.Duration.iloc[y1.argmax()]/0.13525**3+0.00001
  #y = np.sqrt(np.pi)/2/x*special.erf(x)
  #return y*6./y1.max()



##buoni

##fit








 
def IB_Pulse_Thresh(puls,gb,data,Sigma_min):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------

  
  return




def Compare_Beams(puls):
  
  #STUDIARE VALORI E SE METTERE puls.Pulse==0
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=36000)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=3600)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  #count,div = np.histogram(puls.Time[puls.Pulse==0],bins=360)
  #puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=20.])] += 1
  
  sap0 = puls[(puls.SAP==0)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap0['Time_low'] = sap0.Time_c-sap0.dTime
  sap0.sort('Time_low',inplace=True)
  
  sap1 = puls[(puls.SAP==1)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  sap1['Time_low'] = sap1.Time_c-sap1.dTime
  sap1.sort('Time_low',inplace=True)
  
  sap2 = puls[(puls.SAP==2)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
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
  
  puls.Pulse.loc[(puls.SAP==0)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)]=sap0.Pulse
  
  puls.Pulse.loc[(puls.SAP==1)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)]=sap1.Pulse
  
  puls.Pulse.loc[(puls.SAP==2)&(puls.Pulse<=RFI_percent)&(puls.BEAM>12)]=sap2.Pulse
    
  return

