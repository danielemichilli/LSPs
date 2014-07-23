import pandas as pd
import numpy as np

import time

#import pyximport; pyximport.install()
from get_groups import pulse_code

from Parameters import *

def TimeAll(data):
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
    
  data.Time = data.Time + data.DM * delay
   
  return data


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  
  data.sort(['SAP','BEAM','DM'],inplace=True)

  data['Pulse'] = 0  #meglio cosi o ''?
  code = 0
  
  #for ind, event in data.iterrows(): 
  for ind1, beam_group in data.groupby(['SAP','BEAM'],sort=False):
    
    
    time0 = time.clock()
    data['Pulse'] = pulse_code(beam_group.values,beam_group.T.values) 
    print 'Time: ',time.clock()-time0,' s'
    
    print data
    
    data = data[data.Pulse>0]
    
  return data


def Table(data):
    
  gb = data.groupby(['SAP','BEAM','Pulse'],sort=False)  #provare se va piu veloce mettendo il comando in ogni riga

  pulses = data[data.index.isin(gb.Sigma.idxmax())]  #probabilmente esistono modi piu efficienti
  pulses.index = pulses.Pulse
  
  pulses = pulses.ix[:,['SAP','BEAM','DM','Sigma','Time','Downfact']]
  
  pulses['dDM'] = (gb.DM.max().values - gb.DM.min().values) / 2.
  pulses['dTime'] = (gb.Time.max().values - gb.Time.min().values) / 2.
  pulses['Down_Time'] = gb.Down_Time.max().values
  
  #print pulses
  
  return pulses