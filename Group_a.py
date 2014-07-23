import pandas as pd
import numpy as np

import time

from Parameters import *

def TimeAll(data):
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
    
  data.Time = data.Time + data.DM * delay
   
  return data


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  
  data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)

  
  data['Pulse'] = 0  #meglio cosi o ''?
  code = 0
  
  #for ind, event in data.iterrows(): 
  
  
  
  for ind1, beam_group in data.groupby(['SAP','BEAM'],sort=False):
    
    time0 = time.clock()
    
    for ind, event in beam_group.iterrows():
     #provare anche Cython  #provare modificando righe successive
      
#provere con gruppi
#However, the grouping is only an intermediate step; for example, we may want to iterate over each of the patient groups:
#In [135]:
#for patient, group in cdystonia_grouped:
#    print patient
#    print group
#    print




#probabilmente meglio apply 
#def top(df, column, n=5):
#    return df.sort_index(by=column, ascending=False)[:n]
#
#top3segments = segments_merged.groupby('mmsi').apply(top, column='seg_length', n=3)[['names', 'seg_length']]
#top3segments
     
      if event.DM < 40: step = 0.01
      elif event.DM < 140: step = 0.05
      else: step = 0.1
      
      line = data[(data.DM >= event['DM']-STEPS_GROUP*step-ERR_FLOAT) & (data.DM <= event['DM']-step+ERR_FLOAT)]  #probabilmente meglio mettere in range
      cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(event.Down_Time,line.Down_Time)  #forse si puo allargare il range per selezionare rfi 
      line = line[cond]
      line = line[line.DM==line.DM.max()]
      rows = len(line)
      
      if rows > 1:
        erase = line[line.index!=line.Sigma.idxmax()]
        data.drop(erase.index,inplace=True)
        line = line[line.index==line.Sigma.idxmax()]
        rows = len(line)
        
      if rows == 0:
        code += 1
        data.loc[ind,'Pulse'] = code
  
      if rows == 1:
        data.loc[ind,'Pulse'] = line['Pulse'].iloc[0]
  
    print 'Time: ',time.clock()-time0,' s'

  data.drop(data.index[data.Pulse == 0],inplace=True)
  data = data.groupby('Pulse',sort=False).filter(lambda x: len(x) > 3)
  
  return data

  

def Table(data):
   
  
  gb = data.groupby(['SAP','BEAM','Pulse'],sort=False)  #provare se va piu veloce mettendo il comando in ogni riga
  
  pulses = data.reindex(gb.Sigma.idxmax())  #probabilmente esistono modi piu efficienti
  
  pulses.index = pulses.Pulse
  
  pulses = pulses.ix[:,['SAP','BEAM','DM','Sigma','Time','Downfact']]
  
  pulses['DM_min']= gb.DM.min().values
  pulses['Sigma_min'] = gb.Sigma.min().values
  pulses['dDM'] = (gb.DM.max().values - gb.DM.min().values) / 2.
  pulses['dTime'] = (gb.Time.max().values - gb.Time.min().values) / 2.
  pulses['Down_Time'] = gb.Down_Time.max().values
  
  #print pulses
  
  return pulses