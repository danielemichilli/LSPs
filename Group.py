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


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  #------------------------------
  # Assigns each event to a pulse
  #------------------------------
  
  data.sort(['DM'],inplace=True)
  data['Pulse'] = 0
  data.Pulse = data.Pulse.astype(np.int32)
  data_copy = data.ix[:,['DM','Sigma','Time','Duration']]

  C_Funct.Get_Group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Duration.values,data.Pulse.values)
  
  data = data[data.Pulse>0]
  
  return data


#UNIRE le due funzioni
def Table(data):
  #-----------------------------
  # Create a table of the pulses
  #-----------------------------
  
  data.Time = TimeAlign(data.Time,data.DM)
  
  gb = data.groupby('Pulse',sort=False)  #provare se va piu veloce mettendo il comando in ogni riga
  
  #print data[data==9489].groupby('Pulse').Time.value_counts()
  #print data[data==9489].groupby('Pulse').std()

  pulses = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  pulses.index = pulses.Pulse
  pulses = pulses.loc[:,['DM','Sigma','Time','Duration']]
  pulses['Pulse']=1  #METTERE prima
  pulses.Pulse=pulses.Pulse.astype(np.int8)
  
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['DM_min'] = gb.DM.min()
  pulses['Sigma_min'] = gb.Sigma.min()
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  
  pulses['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  pulses.DM_c=pulses.DM_c.astype(np.float32)
  pulses['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  pulses.Time_c=pulses.Time_c.astype(np.float32)
  
  
  #mettere in RFIexcision
  pulses['Sigma_DM_max'] = data.Sigma[gb.DM.idxmax()].values
  pulses['Sigma_DM_min'] = data.Sigma[gb.DM.idxmin()].values
  
  pulses['N_events'] = gb.DM.count()
  pulses.N_events=pulses.N_events.astype(np.int16)
  
  
  pulses = RFIexcision.Group(pulses,gb)
  
  return pulses