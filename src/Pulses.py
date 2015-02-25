import numpy as np
import pandas as pd

import C_Funct
import RFIexcision
from Parameters import *

def Initialize():
  #Creates the tables in memory
  pulses = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse','dDM','dTime','DM_c','Time_c','N_events'])
  
  pulses = pulses.astype(np.float32)
  pulses.SAP = pulses.SAP.astype(np.uint8)
  pulses.BEAM = pulses.BEAM.astype(np.uint8)
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses.N_events = pulses.N_events.astype(np.int16)
  
  return pulses
  

def Generator(events):
  #-------------------------------
  # Create a table with the pulses
  #-------------------------------
  
  pulses = Initialize()
  
  gb = events.groupby('Pulse',sort=False)
  pulses = events.loc[events.index.isin(gb.Sigma.idxmax())]
  pulses.index = pulses.Pulse
  pulses.index.name = None
  pulses = pulses.loc[:,['SAP','BEAM','DM','Sigma','Time','Duration']]
  pulses['Pulse'] = 0
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  pulses['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  pulses.DM_c=pulses.DM_c.astype(np.float32)
  pulses['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  pulses.Time_c=pulses.Time_c.astype(np.float32)
  pulses['N_events'] = gb.DM.count()
  pulses.N_events = pulses.N_events.astype(np.int16)

  #####   MODIFICARE   #########
  
  pulses = pulses[pulses.N_events>4]
  
  # Reduce the RFI and corrects for the time misalignment
  
  #data_idxmax = data.loc[gb.DM.idxmax()]
  #data_idxmax.index = data_idxmax.Pulse
  #data_idxmin = data.loc[gb.DM.idxmin()]
  #data_idxmin.index = data_idxmin.Pulse  
  #Sigma_DM_max = data_idxmax.Sigma
  #Sigma_DM_min = data_idxmin.Sigma
  #Time_DM_max = data_idxmax.Time
  #Time_DM_min = data_idxmin.Time
  Sigma_min = gb.Sigma.min()
  

  RFIexcision.IB_Pulse_Thresh(pulses,gb,events,Sigma_min)
  RFIexcision.Pulse_Thresh(pulses,gb,events,Sigma_min)
  TimeAlign(pulses.Time,pulses.DM)
  TimeAlign(pulses.Time_c,pulses.DM_c)

  RFIexcision.IB_Align_Pulse_Thresh(pulses,gb,events)
  RFIexcision.Align_Pulse_Thresh(pulses,gb,events)
  
  pulses = pulses[pulses.Pulse <= RFI_percent]

  return pulses
  
  
def TimeAlign(Time,DM):
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 0.000241*2  #s-1
  delay1 = np.float32((F_MIN**-2 - F_MAX**-2) / (k))
  
  #bin=(F_MAX-F_MIN)/288  #288 subbands
  f_range = F_MAX-F_MIN
  #delay2 = np.float32(((F_MIN+f_range/2)**-2 - (F_MIN+f_range*145/288)**-2) / k)
  delay2=0
  
  Time += DM * (delay1 + delay2)
  
  return Time