#############################
#
# Group
#
# Written by Daniele Michilli
#
#############################

import pandas as pd
import numpy as np

import C_Funct
import RFIexcision

from Parameters import *


def TimeAlign(Time,DM):
  #----------------------
  # Aligns pulses in time
  #----------------------

  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
  
  Time = Time + DM * delay
  Time = Time.astype(np.float32)
  
  return Time


def Pulses(data):
  #------------------------------
  # Assigns each event to a pulse
  #------------------------------
  
  data.sort(['DM'],inplace=True)
  data['Pulse'] = 0
  data.Pulse = data.Pulse.astype(np.int32)
  data_copy = data.ix[:,['DM','Sigma','Time','Duration']]

  C_Funct.Get_Group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Duration.values,data.Pulse.values)
  
  data = data[data.Pulse>0]
  
  #-----------------------------
  # Create a table of the pulses
  #-----------------------------
  
  data.Time = TimeAlign(data.Time,data.DM)
  
  gb = data.groupby('Pulse',sort=False)
  
  puls = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  puls.index = puls.Pulse
  puls = puls.loc[:,['DM','Sigma','Time','Duration']]
  puls['Pulse']=1
  puls.Pulse=puls.Pulse.astype(np.int8)
  
  puls['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  puls.dDM=puls.dDM.astype(np.float32)
  puls['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  puls.dTime=puls.dTime.astype(np.float32)
  puls['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  puls.DM_c=puls.DM_c.astype(np.float32)
  puls['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  puls.Time_c=puls.Time_c.astype(np.float32)
  
  puls = RFIexcision.Pulse_Thresh(puls,gb,data)
  
  return puls