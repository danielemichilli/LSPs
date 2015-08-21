import pandas as pd
import numpy as np
import multiprocessing as mp

import Pulses


def candidates(pulses,folders):
  #Create candidates lists
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  pool = mp.Pool()
  results = pool.map(Repeated_candidates_beam, [(pulses,sap,beam) for (sap,beam) in folders])
  pool.close()
  pool.join()
  #Pulses.Repeated_candidates(pulses)
  
  pulses.Candidate = pd.concat(results)
  results = 0
  
  #Unify the same repeated candidates in different beams
  pulses.sort('DM',inplace=True)
  for sap in range(3):
    cands_SAP = pulses[(pulses.SAP==sap)&(pulses.Candidate>0)]
    if not cands_SAP.empty:
      diff_DM = np.abs(cands_SAP.DM-cands_SAP.DM.shift())
      diff_DM.iloc[0] = diff_DM.iloc[1]
      diff_DM -= 0.5
      diff_DM = diff_DM.round()
      cands_SAP = cands_SAP.Candidate
      cands_SAP[diff_DM<.1] = np.nan
      cands_SAP.fillna(method='pad',inplace=True)
      pulses.Candidate.loc[cands_SAP.index] = cands_SAP
    
    
  
  
  
  #si puo' aggiungere una funzione per non avere lo stesso unique candidate in beams differenti
  pulses.sort('Sigma',ascending=False,inplace=True)
  cands_unique = pulses[(pulses.Pulse==0)&(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby('N_events')['Sigma'].nlargest(4)
  pulses.Candidate.loc[cands_unique.index.get_level_values(1)] = np.arange(cands_unique.shape[0])
  
  cands = candidates_generator(pulses[pulses.Candidate>=0])
  
  return cands


def candidates_generator(pulses):
  #Prende max N_pulses e mean DM tra tutti i beams, non so se sia l'opzione migliore
  cands = pulses.groupby(['Candidate','BEAM']).agg({'Sigma':np.sum,'SAP':np.size,'DM':np.mean}).groupby(level=0).agg({'Sigma':np.max,'SAP':np.max,'DM':np.mean})
  cands.index.name = None
  cands.rename(columns={'SAP': 'N_pulses'}, inplace=True)




 


def Repeated_candidates_beam((pulses,sap,beam)):
  pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)&(pulses.Pulse<=2)]
  pulses.sort('Sigma',inplace=True)
  pulses.DM = pulses.DM.astype(np.float64).round(2)
  
  n_pulses = pulses.shape[0]
  if n_pulses < 32: min_elements = 2
  elif n_pulses < 115: min_elements = 3
  elif n_pulses < 233: min_elements = 4
  elif n_pulses < 373: min_elements = 5
  elif n_pulses < 526: min_elements = 6
  else: min_elements = 7
  
  span = 0.05
  
  top_count = pulses.groupby('DM')['Sigma'].count()
  top_sum = pulses.groupby('DM')['Sigma'].sum()
  
  top_count = top_count[top_count >= min_elements]
  top_sum = top_sum[top_sum >= min_elements]

  cand = pulses.N_events
  cand[:] = -1
  i = 10

  while not top_sum[top_sum!=0].empty:
    DM = top_sum.argmax()
    Sigma = top_sum.loc[DM-span:DM+span].sum()
    N_pulses = top_count.loc[DM-span:DM+span].count()
    cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    top_count.loc[DM-span:DM+span] = 0
    top_sum.loc[DM-span:DM+span] = 0
    i += 1
    
  return cand





#tolleranza variabile
#def Repeated_candidates_beam((pulses,sap,beam)):  
  #pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)&(pulses.Pulse<=2)]
  
  #n_pulses = pulses.shape[0]
  
  #if n_pulses < 32: min_elements = 2
  #elif n_pulses < 115: min_elements = 3
  #elif n_pulses < 233: min_elements = 4
  #elif n_pulses < 373: min_elements = 5
  #elif n_pulses < 526: min_elements = 6
  #else: min_elements = -1
  
  #if min_elements < 0: span = 0.05
  #else: span = 545.*(10./222./Utilities.c(n_pulses,min_elements)**(1./(min_elements-1)))
  #if span > 0.25: span = 0.25
  
  #pulses.sort('Sigma',inplace=True)
  #pulses.DM = pulses.DM.astype(np.float64).round(2)
  
  #top_count = pulses.groupby('DM')['Sigma'].count()
  #top_sum = pulses.groupby('DM')['Sigma'].sum()
  
  #if min_elements < 0:
    #top_count = top_count.nlargest(10)
    #top_sum = top_sum.loc[top_count.index]
  #else:
    #top_count = top_count[top_count >= min_elements]
    #top_sum = top_sum[top_sum >= min_elements]

  ##cand = pd.DataFrame(columns=['DM','Sigma','N_pulses'])
  #cand = pulses.N_events
  #cand[:] = -1
  #i = 10

  #while not top_sum[top_sum!=0].empty:
    #DM = top_sum.argmax()
    #Sigma = top_sum.loc[DM-span:DM+span].sum()
    #N_pulses = top_count.loc[DM-span:DM+span].count()
    ##cand.loc[i] = ([DM,Sigma,N_pulses])
    ##pulses.Candidate[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    #cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    #top_count.loc[DM-span:DM+span] = 0
    #top_sum.loc[DM-span:DM+span] = 0
    #i += 1
    
  ##cand.index = (cand.index * 10 + sap) * 100 + beam
    
  #return cand

