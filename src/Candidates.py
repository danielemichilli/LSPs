import pandas as pd
import numpy as np
import multiprocessing as mp
import logging

import Pulses
import C_Funct
import Utilities


def candidates(pulses,idL):
  pulses.sort(['Pulse','Sigma'],ascending=[1,0],inplace=True)

  pulses.Candidate[pulses.Pulse==0] = repeated_parallel(pulses[pulses.Pulse==0],0)
  
  if pulses.Candidate.unique().size <=12:
    pulses.Candidate[(pulses.Pulse<=1)&(pulses.Candidate<0)] = repeated_parallel(pulses[(pulses.Pulse<=1)&(pulses.Candidate<0)],1)

  if pulses.Candidate.unique().size <=12:
    pulses.Candidate[pulses.Candidate<0] = repeated_parallel(pulses[pulses.Candidate<0],2)
  
  cands_unique = pulses[(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby(['SAP','BEAM'],sort=False)[['SAP','BEAM']].head(2).astype(np.int16)
  pulses.Candidate.loc[cands_unique.index.get_level_values('idx')] = np.arange(cands_unique.shape[0]) * 1000 + cands_unique.SAP * 100 + cands_unique.BEAM
  
  if not pulses[pulses.Candidate>=0].empty:
    cands = candidates_generator(pulses[pulses.Candidate>=0].copy(),idL)
    cands['main_cand'] = 0
  
    #Unify the same repeated candidates in different beams
    cands.sort(['Rank','Sigma'],ascending=[0,1],inplace=True)
    new_cand = cands.index

    C_Funct.Compare_candidates(cands.DM.astype(np.float32).values,cands.Time.astype(np.float32).values,cands.index.values,cands.main_cand.values)  
  
  else: cands = pd.DataFrame()
  
  return cands


def repeated_parallel(pulses,rank):
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  dirs = [n for n in gb_puls.indices.iterkeys()]
  
  CPUs = mp.cpu_count()
  dirs_range = int(np.ceil(len(dirs)/float(CPUs)))

  pool = mp.Pool()
  dirs_CPU = [dirs[i*dirs_range:(i+1)*dirs_range] for i in range(CPUs) if len(dirs[i*dirs_range:(i+1)*dirs_range]) > 0]
  results = [pool.apply_async(Repeated_candidates_beam, args=(pulses[(pulses.Pulse==0)&pulses.SAP.isin(zip(*dir_CPU)[0])&pulses.BEAM.isin(zip(*dir_CPU)[1])],dir_CPU,rank)) for dir_CPU in dirs_CPU]
  pool.close()
  pool.join()
  return pd.concat([p.get() for p in results])


def Repeated_candidates_beam(pulses,dirs,rank):
  result = pd.Series()
  for (sap,beam) in dirs:
    pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)].copy()
    pulses.DM = 3*(pulses.DM.astype(np.float64)/3).round(2)
    
    span = 0.25
    
    top_count = pulses.groupby('DM')['Sigma'].count()
    top_sum = pulses.groupby('DM')['Sigma'].sum()
    
    top_sum = top_sum[top_count >= 2]
    #top_count = top_count[top_count >= 2]

    cand = pulses.N_events.astype(np.int32)
    cand[:] = -1
    i = 10

    while not top_sum[top_sum!=0].empty:
      DM = top_sum.argmax()
      #Sigma = top_sum.loc[DM-span:DM+span].sum()
      #N_pulses = top_count.loc[DM-span:DM+span].sum()
      if cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)].shape[0] > 1:
        cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i+1) * 10000 + sap * 1000 + beam * 10 + rank
      #top_count.loc[DM-span:DM+span] = 0
      top_sum.loc[DM-span:DM+span] = 0
      i += 1
    
    result = pd.concat((cand,result))
    
  return result


def period(x):
  if x.size<=1: return 0
  else: return Utilities.rrat_period(x)[0]

def period_err(x):
  if x.size<=1: return 0
  else: return Utilities.rrat_period(x)[1]

def candidates_generator(pulses,idL):
  pulses['Period'] = pulses.Time
  pulses['Period_err'] = pulses.Time
  
  cands = pulses.groupby(['Candidate','SAP','BEAM'],as_index=False,sort=False).agg({'Sigma':np.sum,'N_events':np.size,'DM':np.mean,'Time':np.min,'Period':period,'Period_err':period_err,'Pulse':np.max})
  
  cands = cands.astype(np.float32)
  cands[['N_events','Pulse','SAP','BEAM']] = cands[['N_events','Pulse','SAP','BEAM']].astype(np.int16)
  
  cands.index = cands.Candidate.astype(int)
  cands.index.name = 'idx'
  cands.rename(columns={'N_events': 'N_pulses', 'Pulse': 'Rank'}, inplace=True)
  cands = cands.drop('Candidate',axis=1)
  cands.Time[cands.N_pulses>1] = 0
  cands['id'] = idL + '_' + cands.SAP.astype(str) + '_' + cands.BEAM.astype(str) + '_' + cands.index.astype(str)
  
  cands.sort(['Rank','Sigma'],ascending=[1,0],inplace=True)
  best_cands = cands[cands.N_pulses==1].groupby('SAP').head(4)
  best_cands = best_cands.append(cands[cands.N_pulses>1].groupby('BEAM').head(2).groupby('SAP').head(4))
  
  return best_cands
