import numpy as np
import pandas as pd

import C_Funct
import RFIexcision
from Parameters import *
import Utilities


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
  pulses['Candidate'] = -1
  pulses.Candidate = pulses.Candidate.astype(np.int16)
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



  
def Repeated_candidates_beam((pulses,sap,beam)):
  pulses = pulses[(pulses.SAP==sap)&(pulses.BEAM==beam)&(pulses.Pulse<=2)]
  
  n_pulses = pulses.shape[0]
  
  if n_pulses < 32: min_elements = 2
  elif n_pulses < 115: min_elements = 3
  elif n_pulses < 233: min_elements = 4
  elif n_pulses < 373: min_elements = 5
  elif n_pulses < 526: min_elements = 6
  else: min_elements = -1
  
  if min_elements < 0: span = 0.05
  else: span = 545.*(10./222./Utilities.c(n_pulses,min_elements)**(1./(min_elements-1)))
  if span > 0.25: span = 0.25
  
  pulses.sort('Sigma',inplace=True)
  pulses.DM = pulses.DM.astype(np.float64).round(2)
  
  top_count = pulses.groupby('DM')['Sigma'].count()
  top_sum = pulses.groupby('DM')['Sigma'].sum()
  
  if min_elements < 0:
    top_count = top_count.nlargest(10)
    top_sum = top_sum.loc[top_count.index]
  else:
    top_count = top_count[top_count >= min_elements]
    top_sum = top_sum[top_sum >= min_elements]

  #cand = pd.DataFrame(columns=['DM','Sigma','N_pulses'])
  cand = pulses.N_events
  cand[:] = -1
  i = 10

  while not top_sum[top_sum!=0].empty:
    DM = top_sum.argmax()
    Sigma = top_sum.loc[DM-span:DM+span].sum()
    N_pulses = top_count.loc[DM-span:DM+span].count()
    #cand.loc[i] = ([DM,Sigma,N_pulses])
    #pulses.Candidate[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    cand[(pulses.DM>=DM-span)&(pulses.DM<=DM+span)] = (i * np.int16(10) + sap) * np.int16(100) + beam
    top_count.loc[DM-span:DM+span] = 0
    top_sum.loc[DM-span:DM+span] = 0
    i += 1
    
  #cand.index = (cand.index * 10 + sap) * 100 + beam
    
  return cand




