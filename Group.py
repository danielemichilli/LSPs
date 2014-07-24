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

from Parameters import *


def TimeAlign(data):
  #----------------------
  # Aligns pulses in time
  #----------------------

  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
  
  data.Time = data.Time + data.DM * delay
  
  return data


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  #------------------------------
  # Assigns each event to a pulse
  #------------------------------
  
  data.sort(['DM'],inplace=True)
  data['Pulse'] = 0
  data_copy = data.ix[:,['DM','Sigma','Time','Duration']].astype(np.float64)  #provare a vedere i tempi se si toglie questa riga
  
  C_Funct.Get_Group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Duration.values,data.Pulse.values)
  
  data = data[data.Pulse>0]
  
  return data


def Table(data):
  #-----------------------------
  # Create a table of the pulses
  #-----------------------------
  
  gb = data.groupby(['Pulse'],sort=False)  #provare se va piu veloce mettendo il comando in ogni riga

  pulses = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  pulses.index = pulses.Pulse
  pulses = pulses.ix[:,['DM','Sigma','Time','Downfact','Duration','Pulse']]
  
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses['DM_min'] = gb.DM.min()
  pulses['Sigma_min'] = gb.Sigma.min()
  
  return pulses