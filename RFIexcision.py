########################################
#
# Radio Frequency Interferences excision
#
# Written by Daniele Michilli
#
########################################

import numpy as np

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
  #data = data[data.Duration<DURATION_MAX]
  
  #count,div = np.histogram(data.Time,bins=36000)
  #data = data[((data.Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count<=3.*np.median(count)])]  #forse metodi piu efficienti
  
  return data



def IB_Event_Thresh(data):
  #--------------------------------------------------------
  # Applies thresholds to the events in an incoherent beams
  #--------------------------------------------------------

  #Remove low DM
  #data = data[data.DM>DM_MIN]
  
  ##Remove low sigma
  #data = data[data.Sigma>SIGMA_MIN]

  ##Remove high durations
  #data = data[data.Duration<DURATION_MAX]
  
  return data



def Pulse_Thresh(puls,gb,data):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  
  puls.Pulse[puls.DM < 3] = 0
  
  puls.Pulse[puls.dDM > 1.] = 0
  
  puls.Pulse[puls.Sigma < 5.5] = 0
  

  Sigma_DM_max = data.Sigma[gb.DM.idxmax()].values
  Sigma_DM_min = data.Sigma[gb.DM.idxmin()].values
  
  DM_min = gb.DM.min()
  Sigma_min = gb.Sigma.min()
  
  N_events = gb.DM.count()

  puls.Pulse[N_events < 5] = 0
  
  puls.Pulse[gb.Time.apply(np.var) > 0.00005] = 0
  
  puls.Pulse[(puls.DM<=40.) & (puls.dDM/(N_events-1)/0.01 > 0.65)] = 0
  puls.Pulse[(puls.DM>40.) & (puls.DM<=140.) & (puls.dDM/(N_events-1)/0.05 > 0.65)] = 0
  puls.Pulse[(puls.DM>140.) & (puls.dDM/(N_events-1)/0.1 > 0.65)] = 0
  
  puls.Pulse[N_events/puls.Sigma > 4.] = 0
  puls.Pulse[abs(puls.DM-puls.DM_c)/puls.dDM**4 > 1000.] = 0  
  
  puls.Pulse[puls.Sigma/Sigma_min**4 < 0.0043] = 0
  
  puls.Pulse[abs(Sigma_DM_max-Sigma_DM_min) > 2.] = 0
  
  puls.Pulse[abs(1.- (puls.DM-DM_min)/puls.dDM) > .75] = 0
  
  puls.Pulse[puls.Sigma/puls.dDM**4 > 85000.] = 0
  
  return puls



  
def IB_Pulse_Thresh(grouped):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------

  grouped.Pulse[grouped.N_events<4] = 0

  
  #data = data.groupby('Pulse',sort=False).filter(lambda x: len(x) > 3)
  
  #data.drop(data[data.Pulse.isin(grouped[grouped.dDM>3.].index)].index,inplace=True)
  #data.drop(data[data.Pulse.isin(grouped[grouped.dTime>3.*grouped.Duration].index)].index,inplace=True)
  
  #cond1 = (grouped.DM/(grouped.DM_min+grouped.dDM) > 0.095) & (grouped.DM/(grouped.DM_min+grouped.dDM) < 1.005)
  #cond2 = grouped.Sigma/grouped.Sigma_min> 1.1

  #data = data[data.Pulse.isin(grouped[cond1 & cond2].index)]
  
  #AGGIUNGERE analisi dell'andamento del SNR (costante vs piccato)
  
  return grouped


def Compare_IB(coh,incoh_temp,incoh):
  
  CB = coh.ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  IB = incoh_temp.ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  C_Funct.Compare_IB(CB.DM_c.values,CB.dDM.values,CB.Time_c.values,CB.dTime.values,CB.Sigma.values,CB.Pulse.values,\
                     IB.DM_c.values,IB.dDM.values,IB.Time_c.values,IB.dTime.values,IB.Sigma.values,IB.Pulse.values)  
  
  coh.Pulse=CB.Pulse
  incoh.Pulse=IB.Pulse
  
  return


def Compare_Beams(puls):
  
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
  
  Time = puls.Time - puls.DM * delay
  
  count,div = np.histogram(Time,bins=36000)
  puls.Pulse[((Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] = 0  #forse metodi piu efficienti
  
  #print np.count_nonzero(count>=2.*np.median(count[count>0]))
  
  count,div = np.histogram(Time,bins=3600)
  puls.Pulse[((Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] = 0  #forse metodi piu efficienti
  
  #print np.count_nonzero(count>=2.*np.median(count[count>0]))
  
  count,div = np.histogram(Time,bins=360)
  puls.Pulse[((Time-0.01)/(div[1]-div[0])).astype(np.int16).isin(div.argsort()[count>=2.*np.median(count[count>0])])] = 0  #forse metodi piu efficienti
  
  #print np.count_nonzero(count>=2.*np.median(count[count>0]))
  
  Time = 0.
  
  
  sap0 = puls[puls.SAP==0].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  sap1 = puls[puls.SAP==1].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
  
  sap2 = puls[puls.SAP==2].ix[:,['DM_c','dDM','Time_c','dTime','Sigma','Pulse']]
   
        
  C_Funct.Compare_CB(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                     sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values)
  
  C_Funct.Compare_CB(sap0.DM_c.values,sap0.dDM.values,sap0.Time_c.values,sap0.dTime.values,sap0.Sigma.values,sap0.Pulse.values,\
                     sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values)
    
  C_Funct.Compare_CB(sap1.DM_c.values,sap1.dDM.values,sap1.Time_c.values,sap1.dTime.values,sap1.Sigma.values,sap1.Pulse.values,\
                     sap2.DM_c.values,sap2.dDM.values,sap2.Time_c.values,sap2.dTime.values,sap2.Sigma.values,sap2.Pulse.values)
  
  
  puls.Pulse[puls.SAP==0]=sap0.Pulse
  
  puls.Pulse[puls.SAP==1]=sap1.Pulse
  
  puls.Pulse[puls.SAP==2]=sap2.Pulse
    
  return puls

