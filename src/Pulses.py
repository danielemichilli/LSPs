import numpy as np
import pandas as pd

import C_Funct
import RFIexcision
from Parameters import *


def Generator(events):
  #-------------------------------
  # Create a table with the pulses
  #-------------------------------
    
  gb = events.groupby('Pulse',sort=False)
  pulses = events.loc[gb.Sigma.idxmax()]  
  pulses.index = pulses.Pulse
  pulses.index.name = None
  pulses = pulses.loc[:,['SAP','BEAM','DM','Sigma','Time','Duration','Sample']]
  pulses['Pulse'] = 0
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  #pulses['Candidate'] = -1
  #pulses.Candidate = pulses.Candidate.astype(np.int16)
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  
  #pulses['dSample'] = (gb.Sample.max() - gb.Sample.min()) / 2.
  #pulses.dSample = pulses.dSample.astype(np.float32)  
  
  pulses['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  pulses.DM_c=pulses.DM_c.astype(np.float32)
  pulses['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  pulses.Time_c=pulses.Time_c.astype(np.float32)
  pulses['N_events'] = gb.DM.count()
  pulses.N_events = pulses.N_events.astype(np.int16)

  pulses = pulses[pulses.N_events>4]
  
  return pulses
  

def Candidates(pulses):
  cand_num = 0
  candidates = pulses.N_events.copy()
  candidates.iloc[:] = -1
  pulses_clean = pulses[pulses.Pulse==0]
  diff_DM = np.abs(pulses_clean.DM-pulses_clean.DM.shift())
  diff_DM.iloc[0] = diff_DM.iloc[1]
  
  for idx,diff in diff_DM.iteritems():
    if diff > 0.5: cand_num += 1
    candidates.loc[idx] = cand_num
    
  pulses.Candidate = candidates.astype(np.int16)
  
  return
        
