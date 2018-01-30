import pandas as pd
import numpy as np
import logging

import Pulses
import C_Funct
import Utilities


def candidates(pulses,idL):
  pulses.sort(['Sigma','Pulse'],ascending=[0,1],inplace=True)
  
  Repeated_candidates_beam(pulses)
  
  cands_unique = pulses[(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby(['SAP','BEAM'],sort=False)[['SAP','BEAM']].head(5).astype(np.int32)
  pulses.Candidate.loc[cands_unique.index.get_level_values('idx')] = 2 * (np.arange(cands_unique.shape[0]) * 10000 + cands_unique.SAP * 1000 + cands_unique.BEAM)
  #Unique candidates have even ID
  
  if not pulses[pulses.Candidate>=0].empty:
    cands = candidates_generator(pulses[pulses.Candidate>=0].copy(), idL)
    cands['main_cand'] = 0
  
    #Unify the same repeated candidates in different beams
    cands.sort('Sigma',ascending=True,inplace=True)

    C_Funct.Compare_candidates(cands.DM.astype(np.float32).values,cands.Time.astype(np.float32).values,cands.index.values,cands.main_cand.values)
    
    cands.sort(['main_cand', 'Sigma'], ascending=[1,0], inplace=True)
    
  else: cands = pd.DataFrame()
  
  return cands


def Repeated_candidates_beam(pulses):
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  dirs = [n for n in gb_puls.indices.iterkeys()]
  pulses.Candidate[:] = -1

  span = 0.151

  for (sap,beam) in dirs:
    puls_beam = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)]
    
    def group_SNR(DM, pulses):                                                 
      return pulses.Sigma[(pulses.DM >= DM - span) & (pulses.DM <= DM + span)].sum()
    def group_count(DM, pulses):                                                 
      return pulses[(pulses.DM >= DM - span) & (pulses.DM <= DM + span)].shape[0]
    puls_beam['top_SNR'] = puls_beam.apply(lambda x: group_SNR(x.DM, puls_beam), axis=1)
    puls_beam['top_count'] = puls_beam.apply(lambda x: group_count(x.DM, puls_beam), axis=1)
    puls_beam = puls_beam[puls_beam.top_count > 1]
    
    i = 1
    while puls_beam.shape[0] > 0:
      DM = puls_beam.DM[puls_beam.top_SNR.argmax()]
      selected_pulses = puls_beam.Candidate[(puls_beam.DM >= DM - span) & (puls_beam.DM <= DM + span)]
      if selected_pulses.shape[0] > 1:
        pulses.Candidate.loc[selected_pulses.index] = 1 + 2 * (i * 10000 + sap * 1000 + beam)  #Repeated candidates have odd ID
      puls_beam = puls_beam.drop(selected_pulses.index)
      i += 1
    
  return


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
  return cands
