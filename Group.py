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
  
  data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)

  
  data['Pulse'] = ''  #meglio cosi o zeri?
  code = 0
  
  print 'ok'
  
  for ind,event in data.iterrows():  #probabilmente meglio apply
    
    #if event['DM']==500: break
    
    if event.DM < 40: step = 0.01
    elif event.DM < 140: step = 0.05
    else: step = 0.1
    
    #event.Pulse=data.Pulse[ind]
    
    new = (event.Pulse == '')

    
    line = data[(data.Pulse=='') & (data.DM == event.DM+step)]  #(data.DM >= event['DM']+step) & (data.DM <= event['DM']+3.*step)]
    cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(np.multiply(event.Downfact,event.Sampling),np.multiply(line.Downfact,line.Sampling))  #forse meglio selezionare in range
    line = line[cond]
    #line = line[line.DM==line.DM.min()]
    
    
    rows = len(line)
        
    

    if rows > 1:
      erase = line.drop(line[line.Sigma==line.Sigma.max()].index)
      data.drop(erase.index,inplace=True)
      
      line = line[line.Sigma==line.Sigma.max()]
      rows = len(line)
      
      if rows > 1: print 'Attenzione: ',rows
#      
#      line = line[line.Sigma==line.Sigma.max()]
#      rows = len(line)
    if rows == 1:
      if new:
        code += 1
        event['Pulse'] = code
      data.loc[line.index,'Pulse']=event.Pulse
      #data.Pulse[data.index==line.index]=event.Pulse
      
      #if (event.DM < 5) & (event.DM > 4.7): 
      #  print 'event: ',event
      #  print 'line: ',line
      #  print 'data',data.Pulse[data.index==line.index]


    # if event.Pulse == '': data.drop(ind,inplace=True)
    
 
  data.drop(data.index[data.Pulse == ''],inplace=True)
  data = data.groupby('Pulse').filter(lambda x: len(x) > 5)
    
 
      

#    for ind2,event_new in data[data.DM==event['DM']+step].iterrows():
#      if event_new['Pulse'] == '':
#        if abs(event.Time - event_new.Time) <= (2. * (event.Downfact * event.Sampling + event_new.Downfact * event_new.Sampling)):  #Si potrebbe mettere una condizione sulla sigma dei due eventi per RFI
#                                                  #meglio cosi o con gruppi?
#          if event['Pulse'] == '': 
#            code += 1
#            event['Pulse'] = code 
#          event_new['Pulse'] = event['Pulse']
#          
#          break
#    if event['Pulse'] == '': data.drop(ind,inplace=True)
  
  return

      
      
#  pulse_old = pd.DataFrame()
#  for i in np.linspace(0,9.9,0.1):
#    pulse_new = data[(data.DM>i) & (data.DM<=i+0.1)]
#    msk = pulse_new.merge(pulse_old,on='DM',suffixes=['_new','_old'],copy=False,right_index=True)
#    for t in msk.iterrows()