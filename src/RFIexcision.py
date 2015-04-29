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

import C_Funct
from Parameters import *

import time

def Pulse_Thresh(puls,gb,data,Sigma_min):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  
  min_chunk = [5,9,13]
  max_chunk = list(min_chunk)
  max_chunk[0] = np.inf
  max_chunk=np.roll(np.add(max_chunk,-1),-1)
    
  puls_astro = puls[(puls.Pulse < RFI_percent)&(puls.BEAM>12)]

  for i in range(len(min_chunk)):
    puls_chunk = puls_astro[(puls.N_events >= min_chunk[i])&(puls.N_events <= max_chunk[i])].reindex_like(puls)

    puls.Pulse[puls_chunk.Duration > FILTERS[0][i]] += 1
    puls.Pulse[(puls_chunk.DM<=40.48) & (puls_chunk.dDM/(puls.N_events-1)/0.01 > FILTERS[1][i])] += 1
    puls.Pulse[(puls_chunk.DM>40.48) & (puls_chunk.DM<=141.68) & (puls.dDM/(puls.N_events-1)/0.05 > FILTERS[1][i])] += 1
    puls.Pulse[(puls_chunk.DM>141.68) & (puls_chunk.dDM/(puls.N_events-1)/0.1 > FILTERS[1][i])] += 1
    puls.Pulse[abs(puls.DM-puls.DM_c)/puls.dDM > FILTERS[2][i]] += 1  #mettere condizione su N_elements: ignorare se =5 (es. *(N_elements-5))
    puls.Pulse[puls.Sigma/Sigma_min < FILTERS[3][i]] += 1
    puls.Pulse[puls.Sigma/Sigma_min**4 < FILTERS[4][i]] += 1
    
    #Sigma/Sigma_DM_max
    #Sigma/Sigma_DM_min
    
    #puls.Pulse[abs(Sigma_DM_max-Sigma_DM_min) > FILTERS[5][i]] += 1
    #f = data[data.Pulse.isin(puls_chunk.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
    #puls.Pulse[(f>FILTERS[6][i])|(f<FILTERS[7][i])] += 1
    
    #puls.dDM/(puls.dTime+0.0000000001)*(Time_DM_min-Time_DM_max)/(abs(Time_DM_min-Time_DM_max)+0.0000000001)
    
    #STUDIARE nuovamente questi parametri!
    
  
  #puls.Pulse[puls.Duration > 0.0221184] += 3 #BETTER to divide in 3 steps, accordingly to dDMs
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM<=40.5) & (puls.dDM/(N_events-1)/0.01 > 0.667)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM>40.5) & (puls.DM<=141.7) & (puls.dDM/(N_events-1)/0.05 > 0.667)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM>141.7) & (puls.dDM/(N_events-1)/0.1 > 0.667)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (N_events/puls.Sigma > 17.8)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (abs(puls.DM-puls.DM_c)/puls.dDM > 0.714)] += 1  #mettere condizione su N_elements: ignorare se =5 (es. *(N_elements-5))
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.Sigma/Sigma_min < 1.06)] += 1    #*(N_elements-5)/10.
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.Sigma/Sigma_min**4 < 0.00797)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (abs(Sigma_DM_max-Sigma_DM_min) > 1.02)] += 1
  #puls.Pulse[(puls.Pulse<RFI_percent) & (puls.dDM/puls.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max) < 11.3)] += 1  
  #f = data[data.Pulse.isin(puls[puls.Pulse<RFI_percent].index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  #puls.Pulse[(puls.Pulse<RFI_percent) & ((f<-33.6)|(f>-8.55))] += 1




  #min_chunk = [5,9,13]
  #max_chunk = min_chunk
  #max_chunk[0] = np.inf
  #max_chunk=np.roll(np.add(max_chunk,-1),-1)
  
  #puls_astro = puls[(puls.Pulse < RFI_percent)&(puls.BEAM>12)]
  
  #for i in range(len(min_chunk)):
    #puls_chunk = puls_astro[(puls.N_events >= min_chunk[i])&(puls.N_events <= max_chunk[i])]

    #puls.Pulse[gb.Time.apply(np.var) > FILTERS[8][i]] += 1
  
  #studiare se meglio std o senza /N o con altre opzioni
  
  
  
  #Escludere RFI?
  count,div = np.histogram(puls.Time,bins=360)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] += 1
  
  count,div = np.histogram(puls.Time,bins=3600)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=1.*np.median(count[count>0])])] += 1
  
  
  return



  
def IB_Pulse_Thresh(puls,gb,data,Sigma_min):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------
  
  #puls.Pulse[puls.N_events<4] = 10
  
  min_chunk = [4,9]
  max_chunk = list(min_chunk)
  max_chunk[0] = np.inf
  max_chunk=np.roll(np.add(max_chunk,-1),-1)
  
  puls_astro = puls[(puls.Pulse < RFI_percent)&(puls.BEAM==12)]
  
  for i in range(len(min_chunk)):
    puls_chunk = puls_astro[(puls.N_events >= min_chunk[i])&(puls.N_events <= max_chunk[i])].reindex_like(puls)
    
    puls.Pulse[puls_chunk.Duration > FILTERS[0][i]] += 1
    puls.Pulse[(puls_chunk.DM<=40.5) & (puls_chunk.dDM/(puls.N_events-1)/0.01 > FILTERS[1][i])] += 1
    puls.Pulse[(puls_chunk.DM>40.5) & (puls_chunk.DM<=141.7) & (puls.dDM/(puls.N_events-1)/0.05 > FILTERS[1][i])] += 1
    puls.Pulse[(puls_chunk.DM>141.7) & (puls_chunk.dDM/(puls.N_events-1)/0.1 > FILTERS[1][i])] += 1
    puls.Pulse[abs(puls.DM-puls.DM_c)/puls.dDM > FILTERS[2][i]] += 1  #mettere condizione su N_elements: ignorare se =5 (es. *(N_elements-5))
    puls.Pulse[puls.Sigma/Sigma_min < FILTERS[3][i]] += 1



 #min_chunk = [5,9]
  #max_chunk = min_chunk
  #max_chunk[0] = np.inf
  #max_chunk=np.roll(np.add(max_chunk,-1),-1)
  
  #puls_astro = puls[(puls.Pulse < RFI_percent)&(puls.BEAM==12)]
  
  #for i in range(len(min_chunk)):
    #puls_chunk = puls_astro[(puls.N_events >= min_chunk[i])&(puls.N_events <= max_chunk[i])]
    #puls.Pulse[gb.Time.apply(np.var) > FILTERS[6][i]] += 1
  
  #studiare se meglio std o senza /N o con altre opzioni
  
  count,div = np.histogram(puls.Time,bins=360)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] += 1
  
  count,div = np.histogram(puls.Time,bins=3600)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=1.*np.median(count[count>0])])] += 1
  
  
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

