import pandas as pd
import numpy as np

import time


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz

def TimeAll(data):
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
    
  data.Time = data.Time + data.DM * delay
   
  return data


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  
  ERR = 0.001
  STEPS = 2
  
  #data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)

  data['Pulse'] = 0  #meglio cosi o ''?
  code = 0

  #for ind, event in data.iterrows(): 
  for ind1, beam_group in data.groupby(['SAP','BEAM'],sort=False):
    
    time0 = time.clock()
    
    for ind2, DM_group in data.groupby(['DM'],sort=True):
      
      if ind2 < 40: step = 0.01
      elif ind2 < 140: step = 0.05
      else: step = 0.1
      
      line = data[(data.DM >= ind2-STEPS*step-ERR) & (data.DM <= ind2-step+ERR)]
      
      
      
      DM_group.apply(lambda x: group(x,data,code), axis=1)
    
    #beam_group.apply(group,args=(data,code),axis=1) #,raw=True)     
     #provare anche Cython  #provare modificando righe successive
    print '\n',ind1
    print 'Time: ',time.clock()-time0,' s'
  
  #print data

  return data


def group(event,data,code):
      
      cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(event.Down_Time,line.Down_Time)
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
        #data.loc[event.index,'Pulse'] = code
        event.Pulse = code
        
      if rows == 1:
        event.Pulse = line['Pulse'].iloc[0]
        #data.loc[event.index,'Pulse'] = line['Pulse'].iloc[0]


  


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