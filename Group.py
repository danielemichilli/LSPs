import pandas as pd
import numpy as np

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
  
  data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)

  
  data['Pulse'] = 0  #meglio cosi o ''?
  code = 0
  
  #for ind, event in data.iterrows(): 
  for ind1, beam_group in data.groupby(['SAP','BEAM']):
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
      
      line = data[(data.DM >= event['DM']-STEPS*step-ERR) & (data.DM <= event['DM']-step+ERR)]  #probabilmente meglio mettere in range
      cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(np.multiply(np.multiply(event.Downfact,event.Sampling),.5),np.multiply(np.multiply(line.Downfact,line.Sampling),.5))  #forse meglio selezionare in range
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

  data.drop(data.index[data.Pulse == 0],inplace=True)
  data = data.groupby('Pulse').filter(lambda x: len(x) > 3)
  
  print data
  
  return data


def Table(data):
  
  pulses = pd.DataFrame(columns=['SAP','BEAM','Pulse','DM','dDM','dTime','Sigma'])
  
  for ind, beam_group in data.groupby(['SAP','BEAM']):
    
    
    pulses['Sigma'] = data.groupby('Pulse').Sigma.max()
    pulses.index = data[data.index.isin(data.groupby('Pulse').Sigma.idxmax())].Pulses
    pulses['DM'] = data[data.index.isin(data.groupby('Pulse').Sigma.idxmax())].DM
    pulses['DM_max'] = data.groupby('Pulse').DM.max()
    pulses['DM_min'] = data.groupby('Pulse').DM.min()
    pulses['Time'] = data[data.index.isin(data.groupby('Pulse').Sigma.idxmax())].Time
    pulses['Time_max'] = data.groupby('Pulse').Down_Time.max()
    pulses['Time_min'] = data.groupby('Pulse').Down_Time.min()


  #si potrebbe aggiungere downfact*sampling
    
  return pulses