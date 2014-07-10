import pandas as pd
import numpy as np

F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz

def TimeAll(data):
  k = 4.1488078e3  #s
  delay = .5 * k * (F_MIN**-2 - F_MAX**-2)
    
  data.Time = data.Time + data.DM * delay
   
  return


def Pulses(data):  #AGGIUNGERE funzione per beam diversi
  
  ERR = 0.001
  
  data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)

  
  data['Pulse'] = 0  #meglio cosi o zeri?
  code = 0
  
  
  for ind,event in data.iterrows():  #probabilmente meglio apply  #provare anche Cython  #provare DataFrame.itertuples(index=True)  #provare modificando righe successive
    
    if event.DM < 40: step = 0.01
    elif event.DM < 140: step = 0.05
    else: step = 0.1
    
    line = data[(data.DM >= event['DM']-step-ERR) & (data.DM <= event['DM']-step+ERR)]  #probabilmente meglio mettere in range
    #if not line.empty: print 'line:\n',line
    cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(np.multiply(event.Downfact,event.Sampling),np.multiply(line.Downfact,line.Sampling))  #forse meglio selezionare in range
    line = line[cond]
    
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
 
  #data.drop(data.index[data.Pulse == ''],inplace=True)
  #data = data.groupby('Pulse').filter(lambda x: len(x) > 5)
  
  return

      
      
#  pulse_old = pd.DataFrame()
#  for i in np.linspace(0,9.9,0.1):
#    pulse_new = data[(data.DM>i) & (data.DM<=i+0.1)]
#    msk = pulse_new.merge(pulse_old,on='DM',suffixes=['_new','_old'],copy=False,right_index=True)
#    for t in msk.iterrows()