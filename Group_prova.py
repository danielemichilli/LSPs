  for ind,event in data.iterrows():  #probabilmente meglio apply  #provare anche Cython  #provare DataFrame.itertuples(index=True)
    
    if event.DM < 40: step = 0.01
    elif event.DM < 140: step = 0.05
    else: step = 0.1
    
    line = data[(data.Pulse=='') & (data.DM == event.DM-step)]  #(data.DM >= event['DM']+step) & (data.DM <= event['DM']+3.*step)]  #probabilmente meglio guardare righe precendenti piuttosto
    cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(np.multiply(event.Downfact,event.Sampling),np.multiply(line.Downfact,line.Sampling))  #forse meglio selezionare in range
    line = line[cond]
    
    rows = len(line)
    
    if rows > 1:
      erase = line.drop(line[line.index==line.Sigma.idxmax()].index)
      #data.drop(erase.index,inplace=True)
      line = line[line.index==line.Sigma.idxmax()]
      rows = len(line)
    
    if rows == 0:
      code += 1
      data.loc[ind,'Pulse'] = code

    if rows == 1:
      data.loc[ind,'Pulse'] = line.Pulse
 
  data.drop(data.index[data.Pulse == ''],inplace=True)
  data = data.groupby('Pulse').filter(lambda x: len(x) > 5)
  
  return





  for ind,event in data.iterrows():  #probabilmente meglio apply  #provare anche Cython  #provare DataFrame.itertuples(index=True)
      
    if event.DM < 40: step = 0.01
    elif event.DM < 140: step = 0.05
    else: step = 0.1
    
    new = (event.Pulse == '')
    
    line = data[(data.Pulse=='') & (data.DM == event.DM+step)]  #(data.DM >= event['DM']+step) & (data.DM <= event['DM']+3.*step)]  #probabilmente meglio guardare righe precendenti piuttosto
    cond = np.absolute(np.subtract(event.Time,line.Time)) < np.add(np.multiply(event.Downfact,event.Sampling),np.multiply(line.Downfact,line.Sampling))  #forse meglio selezionare in range
    line = line[cond]
    
    rows = len(line)
        
    if rows > 1:
      erase = line.drop(line[line.index==line.Sigma.idxmax()].index)
      #data.drop(erase.index,inplace=True)
      line = line[line.index==line.Sigma.idxmax()]
      rows = len(line)

    if rows == 1:
      if new:
        code += 1
        data.loc[ind,'Pulse'] = code
      data.loc[line.index,'Pulse'] = event.Pulse
      
  data.drop(data.index[data.Pulse == ''],inplace=True)
  data = data.groupby('Pulse').filter(lambda x: len(x) > 5)
  
  return