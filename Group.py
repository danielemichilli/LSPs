import pandas as pd
import numpy as np

from get_groups import pulse_code

from Parameters import *

def TimeAll(data):
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
    
  data.Time = data.Time + data.DM * delay
   
  return data


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  
  data.sort(['SAP','BEAM','DM'],inplace=True)

  data['Pulse'] = 0
  
  #for ind, event in data.iterrows(): 
  for ind1, beam_group in data.groupby(['SAP','BEAM'],sort=False):
    
    data_copy = beam_group.ix[:,['DM','Sigma','Time','Down_Time']].astype(np.float64)
    
    group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Down_Time.values,data.Pulse.values)

    data = data[data.Pulse>0]
    
  return data


def Table(data):
    
  gb = data.groupby(['SAP','BEAM','Pulse'],sort=False)  #provare se va piu veloce mettendo il comando in ogni riga

  pulses = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  pulses.index = pulses.Pulse
  
  pulses = pulses.ix[:,['SAP','BEAM','DM','Sigma','Time','Downfact','Down_Time']]
  
  pulses['dDM'] = (gb.DM.max().values - gb.DM.min().values) / 2.
  pulses['dTime'] = (gb.Time.max().values - gb.Time.min().values) / 2.
  
  #print pulses
  
  return pulses