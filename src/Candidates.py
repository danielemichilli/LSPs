import pandas as pd
import numpy as np
import multiprocessing as mp
import logging

import Pulses
import C_Funct
import Utilities


def candidates(pulses):
  pulses.sort(['Pulse','Sigma'],ascending=[1,0],inplace=True)

  #Create candidates lists
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  pool = mp.Pool()
  results = pool.map(Repeated_candidates_beam, [(pulses[pulses.Pulse==0],n,0) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()
  pulses.Candidate[pulses.Pulse==0] = pd.concat(results)
  results = 0
    
  if pulses.Candidate.unique().size <=12:
    pool = mp.Pool()
    results = pool.map(Repeated_candidates_beam, [(pulses[(pulses.Pulse<=1)&(pulses.Candidate<0)],n,1) for n in gb_puls.indices.iterkeys()])
    pool.close()
    pool.join()
    pulses.Candidate[(pulses.Pulse<=1)&(pulses.Candidate<0)] = pd.concat(results) * 10 + 1
    results = 0

  if pulses.Candidate.unique().size <=12:
    pool = mp.Pool()
    results = pool.map(Repeated_candidates_beam, [(pulses[pulses.Candidate<0],n,2) for n in gb_puls.indices.iterkeys()])
    pool.close()
    pool.join()
    pulses.Candidate[pulses.Candidate<0] = pd.concat(results) * 10 + 1
    results = 0
  
  cands_unique = pulses[(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby(['SAP','BEAM'],sort=False)[['SAP','BEAM']].head(2)
  pulses.Candidate.loc[cands_unique.index.get_level_values('idx')] = (np.arange(cands_unique.shape[0]) * 10 + cands_unique.SAP) * 100 + cands_unique.BEAM
  
  if not pulses[pulses.Candidate>=0].empty:
    cands = candidates_generator(pulses[pulses.Candidate>=0].copy())
    cands['main_cand'] = 0
  
    #Unify the same repeated candidates in different beams
    cands.sort('Sigma',inplace=True)
    new_cand = cands.index

    C_Funct.Compare_candidates(cands.DM.values,cands.Sigma.values,cands.Time.values,cands.N_pulses.values,cands.index.values,cands.main_cand.values)  
  
  else: cands = pd.DataFrame()
  
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
  
  cands = pulses.groupby(['Candidate','SAP','BEAM'],as_index=False,sort=False).agg({'Sigma':np.sum,'N_events':np.size,'DM':np.mean,'Time':np.min,'Period':period,'Period_err':period_err,'Pulse':np.max})
  
  cands = cands.astype(np.float32)
  cands[['N_events','Pulse']] = cands[['N_events','Pulse']].astype(np.int16)
  
  
  cands.index = cands.Candidate.astype(int)
  cands.index.name = 'idx'
  cands.rename(columns={'N_events': 'N_pulses', 'Pulse': 'Rank'}, inplace=True)
  cands = cands.drop('Candidate',axis=1)
  cands.Time[cands.N_pulses>1] = 0
  
  cands.sort(['Rank','Sigma'],ascending=[1,0],inplace=True)
  best_cands = cands[cands.N_pulses==1].groupby('SAP').head(4)
  best_cands = best_cands.append(cands[cands.N_pulses>1].groupby('BEAM').head(2).groupby('SAP').head(4))
  
  return best_cands



 


def Repeated_candidates_beam((pulses,(sap,beam),rank)):
  pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)].copy()
  pulses.DM = 3*(pulses.DM.astype(np.float64)/3).round(2)
  
  span = 0.25
  
  top_count = pulses.groupby('DM')['Sigma'].count()
  top_sum = pulses.groupby('DM')['Sigma'].sum()
  
  top_sum = top_sum[top_count >= 2]
  #top_count = top_count[top_count >= 2]

  cand = pulses.N_events.astype(np.int16)
  cand[:] = -1
  i = 10

  while not top_sum[top_sum!=0].empty:
    DM = top_sum.argmax()
    #Sigma = top_sum.loc[DM-span:DM+span].sum()
    #N_pulses = top_count.loc[DM-span:DM+span].sum()
    if cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)].shape[0] > 1:
      cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = ((i * 10 + sap) * 100 + beam) * 10 + rank
    #top_count.loc[DM-span:DM+span] = 0
    top_sum.loc[DM-span:DM+span] = 0
    i += 1
    
  return cand