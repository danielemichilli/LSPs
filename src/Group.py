#############################
#
# Grouping Functions
#
# Contains the functions to
# operate on the pulses.
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
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  #print Time[DM<=40.]
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 4.1488078e3  #s
  delay = np.float32(.5 * k * (F_MIN**-2 - F_MAX**-2))
  
  Time  += DM * delay
  #Time[DM<=40.] += np.float32(.5) * DM * delay  #DA DUE PULSARS, CONTROLLARE!!!!!!!!!!
  
  return Time


def Pulses(data,sap,beam):
  #---------------------------------------------------------
  # Assigns each event to a pulse and remove isolated events
  #---------------------------------------------------------
  
  #data = data[(data.Time>2874)&(data.Time<2875)&(data.DM>41.5)&(data.DM<44.5)]
  #data = data[(data.Time>2001)&(data.Time<2002)&(data.DM>140)&(data.DM<142)]

  
  data.sort(['DM'],inplace=True)
  data['Pulse'] = 0
  data.Pulse = data.Pulse.astype(np.int32)
  data_copy = data.ix[:,['DM','Sigma','Time','Duration']]

  C_Funct.Get_Group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Duration.values,data.Pulse.values)
  
  data = data[data.Pulse>0]
  
  data.Pulse = (data.Pulse*10+sap)*100+beam  #pulse code deve essere unico: non ho trovato un modo migliore per selezionare eventi quando beam diversi hanno stesso codice

  #------------------------------------------------------
  # Create a table with the characteristics of the pulses
  #------------------------------------------------------
  
  gb = data.groupby('Pulse',sort=False)
  
  puls = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  puls.index = puls.Pulse
  puls.index.name = None
  puls = puls.loc[:,['DM','Sigma','Time','Duration']]
  puls['Pulse']=0
  puls.Pulse=puls.Pulse.astype(np.int8)
  
  puls['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  puls.dDM=puls.dDM.astype(np.float32)
  puls['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  puls.dTime=puls.dTime.astype(np.float32)
  puls['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  puls.DM_c=puls.DM_c.astype(np.float32)
  puls['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  puls.Time_c=puls.Time_c.astype(np.float32)
  puls['N_events'] = gb.DM.count()
  puls.N_events=puls.N_events.astype(np.int16)
  
  # Reduce the RFI and corrects for the time misalignment
  
  if beam == 12: puls = RFIexcision.IB_Pulse_Thresh(puls,gb,data)
  else: puls = RFIexcision.Pulse_Thresh(puls,gb,data)
  
  data.Time = TimeAlign(data.Time,data.DM)
  puls.Time = TimeAlign(puls.Time,puls.DM)
  puls.Time_c = TimeAlign(puls.Time_c,puls.DM_c)
  
  if beam == 12: puls = RFIexcision.IB_Align_Pulse_Thresh(puls,gb,data)
  else: puls = RFIexcision.Align_Pulse_Thresh(puls,gb,data)
  
  return data, puls
