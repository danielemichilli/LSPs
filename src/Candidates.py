import pandas as pd
import numpy as np
import multiprocessing as mp
import logging

import Pulses
import C_Funct
import Utilities


def candidates(pulses):
  pulses.sort('Sigma',ascending=False,inplace=True)

  #Create candidates lists
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  pool = mp.Pool()
  results = pool.map(Repeated_candidates_beam, [(pulses[(pulses.Pulse==0)&(pulses.Sigma>=6.5)],n) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()
  pulses.Candidate[(pulses.Pulse==0)&(pulses.Sigma>=6.5)] = pd.concat(results)

  pool = mp.Pool()
  results = pool.map(Repeated_candidates_beam, [(pulses[(pulses.Pulse==0)&(pulses.Candidate<0)],n) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()
  pulses.Candidate[(pulses.Pulse==0)&(pulses.Candidate<0)] = pd.concat(results)

  #pool = mp.Pool()
  #results = pool.map(Repeated_candidates_beam, [(pulses[(pulses.Pulse<=2)&(pulses.Sigma>=6.5)&(pulses.Candidate<0)],n) for n in gb_puls.indices.iterkeys()])
  #pool.close()
  #pool.join()
  #pulses.Candidate[(pulses.Pulse<=2)&(pulses.Sigma>=6.5)&(pulses.Candidate<0)] = pd.concat(results)
  #results = 0
  
  cands_unique = pulses[(pulses.Pulse==0)&(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby(['SAP','BEAM'],sort=False)[['SAP','BEAM']].head(3)
  pulses.Candidate.loc[cands_unique.index.get_level_values('idx')] = (np.arange(cands_unique.shape[0]) * 10 + cands_unique.SAP) * 100 + cands_unique.BEAM

  
  cands = candidates_generator(pulses[pulses.Candidate>=0].copy())
  cands['main_cand'] = 0
  
  #Unify the same repeated candidates in different beams
  cands.sort('Sigma',inplace=True)
  new_cand = cands.index
  C_Funct.Compare_candidates(cands.DM.values,cands.Sigma.values,cands.Time.values,cands.N_pulses.values,cands.index.values,cands.main_cand.values)  
  
  return cands



def period(x):
  if x.size<=1: return 0
  else: return Utilities.rrat_period(x)[0]

def period_err(x):
  if x.size<=1: return 0
  else: return Utilities.rrat_period(x)[1]

def candidates_generator(pulses):
  pulses['Period'] = pulses.Time
  pulses['Period_err'] = pulses.Time
  
  cands = pulses.groupby(['Candidate','SAP','BEAM'],as_index=False).agg({'Sigma':np.sum,'N_events':np.size,'DM':np.mean,'Time':np.min,'Period':period,'Period_err':period_err})

  cands.index = cands.Candidate.astype(int)
  cands.index.name = 'idx'
  cands.rename(columns={'N_events': 'N_pulses'}, inplace=True)
  cands = cands.drop('Candidate',axis=1)
  cands.Time[cands.N_pulses>1] = 0
  return cands



 


def Repeated_candidates_beam((pulses,(sap,beam))):
  pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)]
  pulses.DM = pulses.DM.astype(np.float64).round(2)
  
  n_pulses = pulses.shape[0]
  #Values calculated imposing Utilities.p(n,k) < 10./222 (10 events per observation)
  if n_pulses < 23: min_elements = 2
  elif n_pulses < 202: min_elements = 3
  elif n_pulses < 649: min_elements = 4
  elif n_pulses < 1369: min_elements = 5
  else: min_elements = 6
  
  span = 0.5
  
  top_count = pulses.groupby('DM')['Sigma'].count()
  top_sum = pulses.groupby('DM')['Sigma'].sum()
  
  top_sum = top_sum[top_count >= min_elements]
  top_count = top_count[top_count >= min_elements]

  cand = pulses.N_events
  cand[:] = -1
  i = 10

  while not top_sum[top_sum!=0].empty:
    DM = top_sum.argmax()
    Sigma = top_sum.loc[DM-span:DM+span].sum()
    N_pulses = top_count.loc[DM-span:DM+span].sum()
    cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    top_count.loc[DM-span:DM+span] = 0
    top_sum.loc[DM-span:DM+span] = 0
    i += 1
    
  return cand