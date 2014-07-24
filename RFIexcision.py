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
  data = data[data.Duration<DURATION_MAX]
  
  return data



def IB_Event_Thresh(data):
  #--------------------------------------------------------
  # Applies thresholds to the events in an incoherent beams
  #--------------------------------------------------------

  #Remove low DM
  data = data[data.DM>DM_MIN]
  
  #Remove low sigma
  data = data[data.Sigma>SIGMA_MIN]

  #Remove high durations
  data = data[data.Duration<DURATION_MAX]
  
  return data



def Pulse_Thresh(data,grouped):
  #-----------------------------------------------------
  # Applies thresholds to the pulses in a coherent beams
  #-----------------------------------------------------
  
  
  data = data[data.Pulse>0]
  
  data = data.groupby('Pulse',sort=False).filter(lambda x: len(x) > 3)
    
  data.drop(data[data.Pulse.isin(grouped[grouped.dDM>3.].index)].index,inplace=True)
  data.drop(data[data.Pulse.isin(grouped[grouped.dTime>3.*grouped.Duration].index)].index,inplace=True)
  
  cond1 = (grouped.DM/(grouped.DM_min+grouped.dDM) > 0.095) & (grouped.DM/(grouped.DM_min+grouped.dDM) < 1.005)
  cond2 = grouped.Sigma/grouped.Sigma_min> 1.1

  data = data[data.Pulse.isin(grouped[cond1 & cond2].index)]
  
  #AGGIUNGERE analisi dell'andamento del SNR (costante vs piccato)
  
  return data



  
def IB_Pulse_Thresh(data,grouped):
  #--------------------------------------------------------
  # Applies thresholds to the pulses in an incoherent beams
  #--------------------------------------------------------


  data = data.groupby('Pulse',sort=False).filter(lambda x: len(x) > 3)
  
  data.drop(data[data.Pulse.isin(grouped[grouped.dDM>3.].index)].index,inplace=True)
  data.drop(data[data.Pulse.isin(grouped[grouped.dTime>3.*grouped.Duration].index)].index,inplace=True)
  
  cond1 = (grouped.DM/(grouped.DM_min+grouped.dDM) > 0.095) & (grouped.DM/(grouped.DM_min+grouped.dDM) < 1.005)
  cond2 = grouped.Sigma/grouped.Sigma_min> 1.1

  data = data[data.Pulse.isin(grouped[cond1 & cond2].index)]
  
  #AGGIUNGERE analisi dell'andamento del SNR (costante vs piccato)
  
  return data


def Compare_IB(data,pulses):
  
  a=0
  
  return data


def Compare_Beams(data,puls):
  

  #data = data[data.Pulse.isin(grouped.index)]
  
  #data.drop(data[data.Pulse.isin(grouped[cond]).index].index,inplace=True)
  
  puls_copy = puls
  
  #provare mettendo gb0=puls[puls.SAP==0].groupby('BEAM',sort=False), gb1=... e veder se e' piu' veloce
  
  for ind0, beam0 in puls[puls.SAP==0].groupby('BEAM',sort=False):
    
    beam0_ind = beam0.Pulse.astype(np.int64)
    beam0 = beam0.ix[:,['DM','dDM','Time','Duration','Sigma']].astype([np.float64)
    
    
    for ind1, beam1 in puls[puls.SAP==1].groupby('BEAM',sort=False):
      
      beam1_ind = beam1.Pulse.astype(np.int64)
      beam1 = beam1.ix[:,['DM','dDM','Time','Duration','Sigma']].astype(np.float64)
      
      msk = Compare(beam0.DM.values,beam0.dDM.values,beam0.Time.values,beam0.Duration.values,beam0.Sigma.values,beam0.Pulse.values,\
        beam1.DM.values,beam1.dDM.values,beam1.Time.values,beam1.Duration.values,beam1.Sigma.values,beam1.Pulse.values)
      
    for ind2, beam2 in puls[puls.SAP==2].groupby('BEAM',sort=False):
      
      C_Funct.Compare(beam0,beam2)
      
      
  for ind1, beam1 in puls[puls.SAP==1].groupby('BEAM',sort=False):
    
    for ind2, beam2 in puls[puls.SAP==2].groupby('BEAM',sort=False):
    
      C_Funct.Compare(beam1,beam2)
      

  data.drop(data.index[data.Pulse.isin(puls.index[puls.Pulse==0])],inplace=True)
  puls = puls[puls.Pulse>0]
  #oppure data = data[data.Pulse.isin(pulses.index)]
      
  
  #for ind1, sap_group in group_temp[group_temp.SAP.isin([0,1])].groupby('SAP'):
  
    #for ind2, beam_group in sap_group.groupby('BEAM'):
      
      #for ind3, beam_new in group_temp[group_temp.SAP==ind1+1].groupby('BEAM'):
        
        #beam_group.apply(compare,args=(beam_new,grouped),axis=1) #,raw=True)
        
      #for ind3, beam_new in group_temp[group_temp.SAP==ind1+2].groupby('BEAM'):
        
        #beam_group.apply(compare,args=(beam_new,grouped),axis=1) #,raw=True)    
        
  #data = data[data.Pulse]

  return data

