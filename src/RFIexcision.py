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


def Event_Thresh(data):
  #-----------------------------------------------------
  # Applies thresholds to the events in a coherent beams
  #-----------------------------------------------------

  #Remove low DM
  data = data[data.DM>DM_MIN]
  
  #Remove low sigma
  data = data[data.Sigma>SIGMA_MIN]

  #Remove high durations
  data = data[data.Duration<DURATION_MAX]  #BETTER to divide in 3 steps, accordingly to dDMs
  
  #count,div = np.histogram(data.Time,bins=36000)
  #data = data[((data.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count<=3.*np.median(count)])]  #forse metodi piu efficienti
  
  return data



def IB_Event_Thresh(data):
  #--------------------------------------------------------
  # Applies thresholds to the events in an incoherent beams
  #--------------------------------------------------------

  #Remove low DM
  data = data[data.DM>DM_MIN]
  
  ##Remove low sigma
  #data = data[data.Sigma>SIGMA_MIN]

  ##Remove high durations
  #data = data[data.Duration<DURATION_MAX]
  
  return data



def Pulse_Thresh(puls,gb,data):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  
  #puls.Pulse[puls.dDM > 1.] += 3
  
  Sigma_DM_max = data.Sigma[gb.DM.idxmax()].values
  Sigma_DM_min = data.Sigma[gb.DM.idxmin()].values
  
  Time_DM_max = data.Time[gb.DM.idxmax()].values
  Time_DM_min = data.Time[gb.DM.idxmin()].values
  
  DM_min = gb.DM.min()
  Sigma_min = gb.Sigma.min()
  
  N_events = gb.DM.count()
  
  puls.Pulse[N_events < 5] = 10
  
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM<=40.5) & (puls.dDM/(N_events-1)/0.01 > 0.667)] += 1
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM>40.5) & (puls.DM<=141.7) & (puls.dDM/(N_events-1)/0.05 > 0.667)] += 1
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.DM>141.7) & (puls.dDM/(N_events-1)/0.1 > 0.667)] += 1
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (N_events/puls.Sigma > 17.8)] += 1
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (abs(puls.DM-puls.DM_c)/puls.dDM > 0.714)] += 1  #mettere condizione su N_elements: ignorare se =5 (es. *(N_elements-5))
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.Sigma/Sigma_min < 1.06)] += 1    #*(N_elements-5)/10.
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.Sigma/Sigma_min**4 < 0.00797)] += 1
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (abs(Sigma_DM_max-Sigma_DM_min) > 1.02)] += 1
    
  puls.Pulse[(puls.Pulse<RFI_percent) & (puls.dDM/puls.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max) < 11.3)] += 1  
  
  
  f = data[data.Pulse.isin(puls[puls.Pulse<RFI_percent].index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  puls.Pulse[(puls.Pulse<RFI_percent) & ((f<-33.6)|(f>-8.55))] += 1


  return puls

  
def IB_Pulse_Thresh(puls,gb,data):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------

  N_events = gb.DM.count()

  puls.Pulse[N_events<4] = 10

  
  return puls


def Align_Pulse_Thresh(puls,gb,data):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  
  puls.Pulse[(puls.Pulse<RFI_percent) & (gb.Time.apply(np.var) > 0.000107)] += 1
  
  #studiare se meglio std o senza /N o con altre opzioni
  
  
  count,div = np.histogram(puls.Time,bins=3600)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=1.*np.median(count[count>0])])] += 1  #forse metodi piu efficienti
  
  count,div = np.histogram(puls.Time,bins=360)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] += 1  #forse metodi piu efficienti
  

  
  return puls



def IB_Align_Pulse_Thresh(puls,gb,data):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------

 
  return puls



def Compare_IB(coh,incoh):
  
  CB = coh.ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  IB = incoh.ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  C_Funct.Compare(CB.DM_c.values,CB.dDM.values,CB.Time_c.values,CB.dTime.values,CB.Sigma.values,CB.Pulse.values,\
                  IB.DM_c.values,IB.dDM.values,IB.Time_c.values,IB.dTime.values,IB.Sigma.values,IB.Pulse.values,np.int8(0))  
  
  coh.Pulse=CB.Pulse
  incoh.Pulse=IB.Pulse
  
  return


def Compare_Beams(puls):
  
  #STUDIARE VALORI E SE METTERE puls.Pulse==0
  count,div = np.histogram(puls.Time[puls.Pulse==0],bins=36000)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=1.*np.median(count[count>0])])] += 1
  
  count,div = np.histogram(puls.Time[puls.Pulse==0],bins=3600)
  puls.Pulse[((puls.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] += 1
  

  sap0 = puls[puls.SAP==0].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  sap1 = puls[puls.SAP==1].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  sap2 = puls[puls.SAP==2].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  
  logging.info('Comparison is starting')
  
        
  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,np.int8(1))
  
  logging.info('1/3 completed')
  
  C_Funct.Compare(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  logging.info('2/3 completed')
  
  C_Funct.Compare(sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,\
                  sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values,np.int8(1))
  
  
  puls.Pulse[puls.SAP==0]=sap0.Pulse
  
  puls.Pulse[puls.SAP==1]=sap1.Pulse
  
  puls.Pulse[puls.SAP==2]=sap2.Pulse
    
  return puls




def best_pulses(puls,data):
  
  best = puls  #[puls.Duration<0.0221184]  No duration filter for FRB!!
  data = data[data.Pulse.isin(best.index)]
  gb = data.groupby('Pulse',sort=False)
  
  
  best = best[gb.Time.apply(np.var) < 7.97e-5]
  data = data[data.Pulse.isin(best.index)]
  
  k = 4.1488078e3  #s
  delay = np.float32(.5 * k * (F_MIN**-2 - F_MAX**-2))
  best.Time -= best.DM * delay
  data.Time -= data.DM * delay

  gb = data.groupby('Pulse',sort=False)
  N_events = gb.DM.count()
  
  best0 = best[(best.DM<=40.5) & (best.dDM/(N_events-1)/0.01 < 0.6)]
  best1 = best[(best.DM>40.5) & (best.DM<=141.7) & (best.dDM/(N_events-1)/0.05 < 0.6)]
  best2 = best[(best.DM>141.7) & (best.dDM/(N_events-1)/0.1 < 0.6)]
  
  best = pd.concat([best0,best1,best2])
  
  best0=0
  best1=0
  best2=0
  
  
  data = data[data.Pulse.isin(best.index)]
  gb = data.groupby('Pulse',sort=False)
  N_events = gb.DM.count()
  
  best = best[N_events/best.Sigma < 16.7]
  
  best = best[abs(best.DM-best.DM_c)/best.dDM < 0.625] #mettere condizione su N_elements: ignorare se =5 (es. *(N_elements-5))

  data = data[data.Pulse.isin(best.index)]
  gb = data.groupby('Pulse',sort=False)
  Sigma_min = gb.Sigma.min()

  best = best[best.Sigma/Sigma_min > 1.09]  #*(N_elements-5)/10.

  #data = data[data.Pulse.isin(best.index)]
  #gb = data.groupby('Pulse',sort=False)
  #Sigma_min = gb.Sigma.min()

  #best = best[best.Sigma/Sigma_min**4 > 0.00825]

  data = data[data.Pulse.isin(best.index)]
  gb = data.groupby('Pulse',sort=False)
  Sigma_DM_max = data.Sigma[gb.DM.idxmax()].values
  Sigma_DM_min = data.Sigma[gb.DM.idxmin()].values
  
  best = best[abs(Sigma_DM_max-Sigma_DM_min) < 0.73]

  data = data[data.Pulse.isin(best.index)]
  gb = data.groupby('Pulse',sort=False)
  Time_DM_max = data.Time[gb.DM.idxmax()].values
  Time_DM_min = data.Time[gb.DM.idxmin()].values

  best = best[best.dDM/best.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max) > 13.4]  

  if not best.empty:
    data = data[data.Pulse.isin(best.index)]
    f = data.loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
    best = best[(f>-25.5) & (f<-11.8)]
  
  best.Time  += best.DM * delay
  
  best.sort('Sigma',ascending=False,inplace=True)
  
  return best.groupby(['SAP','BEAM'],sort=False).head(10)
  