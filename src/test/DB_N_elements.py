import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mpl

#FARE: vedere percentuali assolute e relative con rfi

store = pd.HDFStore('N_elements.h5','w')

par_names=['Duration','dDM:(N_events-1):step','N_events:Sigma','abs(DM-DM_c):dDM','abs(DM-DM_c):dDM**4','Sigma:Sigma_min','Sigma:Sigma_min**4','abs(Sigma_DM_max-Sigma_DM_min)','dDM:dTime','f (>)','f (<)','var']

obs_name = []

for path, subdirs, files in os.walk(os.getcwd()):
  for f in files:
    if f.endswith('.hdf5'):
      
      idL = os.path.split(os.path.split(path)[0])[1]

      if idL == 'L196702' : DM = 17.27
      elif idL == 'L196708' : DM = 31.95
      elif idL == 'L202425' : DM = 41.48
      elif idL == 'L203651' : DM = 25.18
      elif idL == 'L204712' : DM = 12.50
      elif idL == 'L204720' : DM = 26.72
      elif idL == 'L206791' : DM = 22.55
      elif idL == 'L214462' : DM = 30.54
      elif idL == 'L232327' : DM = 31.04
      elif idL == 'L94172'  : DM = 34.96
      elif idL == 'L98710'  : DM = 10.67
      elif idL == 'L196709' : DM = 36.36
      elif idL == 'L196710' : DM = 43.50
      elif idL == 'L202429' : DM = 20.88
      elif idL == 'L202433' : DM = 26.30
      elif idL == 'L215819' : DM = 43.49

      if DM < 40.6: step=0.01
      elif DM < 141.7: step=0.05
      else: step=0.1

      DM_hi = DM + step + 0.001
      DM_lo = DM - step - 0.001

      down_DM_hi = DM_hi - 2.
      down_DM_lo = DM_lo - 2.
      up_DM_hi = DM_hi + 2.
      up_DM_lo = DM_lo + 2.


      puls = pd.read_hdf(path+'/SinlgePulses.hdf5',idL+'_pulses')
      puls = puls[(puls.DM>DM_lo)&(puls.DM<DM_hi)]

      num = len(puls)

      data = pd.read_hdf(path+'/SinlgePulses.hdf5',idL)
      data = data[data.Pulse.isin(puls.index)]
      gb = data.groupby('Pulse',sort=False)

      N_events = gb.DM.count()
      puls_all=puls[N_events>=5]
      data_all = data[data.Pulse.isin(puls_all.index)]
      gb = data_all.groupby('Pulse',sort=False)
      N_events = gb.DM.count()

      data60 = pd.DataFrame()
      data75 = pd.DataFrame()
      data95 = pd.DataFrame()
      data99 = pd.DataFrame()
      
      puls_sum = pd.DataFrame()
      
      for n in range(5,40):
        line60 = np.empty(12)
        line75 = np.empty(12)
        line95 = np.empty(12)
        line99 = np.empty(12)

        puls = puls_all[N_events==n]
     
        if not puls_sum.empty:
          puls = puls.append(puls_sum)
        
        num = len(puls)
        
        if num < 10:
          #puls_sum = puls
          continue
        
        puls_sum = pd.DataFrame()
        
        data = data_all[data_all.Pulse.isin(puls.index)]
        gb = data.groupby('Pulse',sort=False)
        Sigma_DM_max = data.Sigma[gb.DM.idxmax()].values
        Sigma_DM_min = data.Sigma[gb.DM.idxmin()].values
        Time_DM_max = data.Time[gb.DM.idxmax()].values
        Time_DM_min = data.Time[gb.DM.idxmin()].values
        DM_min = gb.DM.min()
        Sigma_min = gb.Sigma.min()
        
        puls['parameter'] = puls.Duration
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[0] = puls.parameter.iloc[ int(.6*num) ]
        line75[0] = puls.parameter.iloc[ int(.75*num) ]
        line95[0] = puls.parameter.iloc[ int(.95*num) ]
        line99[0] = puls.parameter.iloc[ int(.99*num) ]
        
          
        puls['parameter'] = puls.dDM/(n-1)/step
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[1] = puls.parameter.iloc[ int(.6*num) ]
        line75[1] = puls.parameter.iloc[ int(.75*num) ]
        line95[1] = puls.parameter.iloc[ int(.95*num) ]
        line99[1] = puls.parameter.iloc[ int(.99*num) ]
        

        puls['parameter'] = n/puls.Sigma
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[2] = puls.parameter.iloc[ int(.6*num) ]
        line75[2] = puls.parameter.iloc[ int(.75*num) ]
        line95[2] = puls.parameter.iloc[ int(.95*num) ]
        line99[2] = puls.parameter.iloc[ int(.99*num) ]
        
        
        puls['parameter'] = abs(puls.DM-puls.DM_c)/puls.dDM
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[3] = puls.parameter.iloc[ int(.6*num) ]
        line75[3] = puls.parameter.iloc[ int(.75*num) ]
        line95[3] = puls.parameter.iloc[ int(.95*num) ]
        line99[3] = puls.parameter.iloc[ int(.99*num) ]
        
        
        puls['parameter'] = abs(puls.DM-puls.DM_c)/puls.dDM**4
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[4] = puls.parameter.iloc[ int(.6*num) ]
        line75[4] = puls.parameter.iloc[ int(.75*num) ]
        line95[4] = puls.parameter.iloc[ int(.95*num) ]
        line99[4] = puls.parameter.iloc[ int(.99*num) ]
        
        
        puls['parameter'] = puls.Sigma/Sigma_min
        sign = '<'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[5] = puls.parameter.iloc[ int(.6*num) ]
        line75[5] = puls.parameter.iloc[ int(.75*num) ]
        line95[5] = puls.parameter.iloc[ int(.95*num) ]
        line99[5] = puls.parameter.iloc[ int(.99*num) ]
        
        
        puls['parameter'] = puls.Sigma/Sigma_min**4
        sign = '<'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[6] = puls.parameter.iloc[ int(.6*num) ]
        line75[6] = puls.parameter.iloc[ int(.75*num) ]
        line95[6] = puls.parameter.iloc[ int(.95*num) ]
        line99[6] = puls.parameter.iloc[ int(.99*num) ]
        

        puls['parameter'] = abs(Sigma_DM_max-Sigma_DM_min)
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[7] = puls.parameter.iloc[ int(.6*num) ]
        line75[7] = puls.parameter.iloc[ int(.75*num) ]
        line95[7] = puls.parameter.iloc[ int(.95*num) ]
        line99[7] = puls.parameter.iloc[ int(.99*num) ]
        
        
        puls['parameter'] = puls.dDM/(puls.dTime+0.0000000001)*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max+0.0000000001)
        sign = '<'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[8] = puls.parameter.iloc[ int(.6*num) ]
        line75[8] = puls.parameter.iloc[ int(.75*num) ]
        line95[8] = puls.parameter.iloc[ int(.95*num) ]
        line99[8] = puls.parameter.iloc[ int(.99*num) ]
        

        f = data[data.Pulse.isin(puls.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
        if f.empty: f=float('NaN')
        
        puls['parameter'] = f
        
        puls = puls.sort('parameter',ascending=True)

        line60[9] = puls.parameter.iloc[ int(.6*num) ]
        line75[9] = puls.parameter.iloc[ int(.75*num) ]
        line95[9] = puls.parameter.iloc[ int(.95*num) ]
        line99[9] = puls.parameter.iloc[ int(.99*num) ]
        
        puls = puls.sort('parameter',ascending=False)

        line60[10] = puls.parameter.iloc[ int(.6*num) ]
        line75[10] = puls.parameter.iloc[ int(.75*num) ]
        line95[10] = puls.parameter.iloc[ int(.95*num) ]
        line99[10] = puls.parameter.iloc[ int(.99*num) ]
        

        F_MIN = 119.43  #MHz
        F_MAX = F_MIN + 31.64  #MHz
        k = 4.1488078e3  #s
        delay = np.float32(.6 * k * (F_MIN**-2 - F_MAX**-2))
        data.Time += data.DM * delay
        puls.Time += puls.DM * delay
        gb = data.groupby('Pulse',sort=False)
        
        puls['parameter'] = gb.Time.apply(np.var)
        sign = '>'

        if sign =='<': 
          puls = puls.sort('parameter',ascending=False)
        elif sign =='>': 
          puls = puls.sort('parameter',ascending=True)

        line60[11] = puls.parameter.iloc[ int(.6*num) ]
        line75[11] = puls.parameter.iloc[ int(.75*num) ]
        line95[11] = puls.parameter.iloc[ int(.95*num) ]
        line99[11] = puls.parameter.iloc[ int(.99*num) ]
        
        line60 = pd.DataFrame(line60).T
        line75 = pd.DataFrame(line75).T
        line95 = pd.DataFrame(line95).T
        line99 = pd.DataFrame(line99).T
        
        line60.index=[n,]
        line75.index=[n,]
        line95.index=[n,]
        line99.index=[n,]
        
        data60=data60.append(line60)
        data75=data75.append(line75)
        data95=data95.append(line95)
        data99=data99.append(line99)
        
      if (not data60.empty) & (not data75.empty) & (not data95.empty) & (not data99.empty) :
        
        obs_name.append(idL)

        data60.columns=[i+'_60' for i in par_names]
        data75.columns=[i+'_75' for i in par_names]
        data95.columns=[i+'_95' for i in par_names]
        data99.columns=[i+'_99' for i in par_names]

        #data60.index += 5
        #data75.index += 5
        #data95.index += 5
        #data99.index += 5
            
        #data60 = data60/data60.loc[data60.abs().idxmax()]
        #data75 = data75/data75.loc[data75.abs().idxmax()]
        #data95 = data95/data95.loc[data95.abs().idxmax()]
        #data99 = data99/data99.loc[data99.abs().idxmax()]
        
        data_store = pd.concat([data60,data75,data95,data99],axis=1)
        store['{}'.format(idL)] = data_store
        

for par in par_names:
  par60 = par+'_60'
  par75 = par+'_75'
  par95 = par+'_95'
  par99 = par+'_99'

  data_plot_60 = pd.DataFrame()
  data_plot_75 = pd.DataFrame()
  data_plot_95 = pd.DataFrame()
  data_plot_99 = pd.DataFrame()
    
  for obs in store.items(): 
    data_plot_60 = pd.concat([data_plot_60,store['{}'.format(obs[0])][par60]],axis=1)
    data_plot_75 = pd.concat([data_plot_75,store['{}'.format(obs[0])][par75]],axis=1)
    data_plot_95 = pd.concat([data_plot_95,store['{}'.format(obs[0])][par95]],axis=1)
    data_plot_99 = pd.concat([data_plot_99,store['{}'.format(obs[0])][par99]],axis=1)

  
    #data_plot_60.fillna(method='bfill',inplace=True)
    #data_plot_75.fillna(method='bfill',inplace=True)
    #data_plot_95.fillna(method='bfill',inplace=True)
    #data_plot_99.fillna(method='bfill',inplace=True)
  

  data_plot_60.columns=obs_name
  data_plot_75.columns=obs_name
  data_plot_95.columns=obs_name
  data_plot_99.columns=obs_name

  mpl.rc('font',size=5)
  
  ax = data_plot_60.plot(title=par60,marker='+')
  ax.yaxis.grid(True, which='minor') 
  ax.set_xlabel("Number of elements")
  ax.set_ylabel("Parameter value")
  plt.minorticks_on()
  plt.grid(which='minor')
  plt.savefig('N_elements/{}.png'.format(par60),format='png',bbox_inches='tight',dpi=200)
  ax = data_plot_75.plot(title=par75,marker='+')
  ax.yaxis.grid(True, which='minor') 
  ax.set_xlabel("Number of elements")
  ax.set_ylabel("Parameter value")
  plt.minorticks_on()
  plt.grid(which='minor')
  plt.savefig('N_elements/{}.png'.format(par75),format='png',bbox_inches='tight',dpi=200)
  ax = data_plot_95.plot(title=par95,marker='+')
  ax.yaxis.grid(True, which='minor') 
  ax.set_xlabel("Number of elements")
  ax.set_ylabel("Parameter value")
  plt.minorticks_on()
  plt.grid(which='minor')
  plt.savefig('N_elements/{}.png'.format(par95),format='png',bbox_inches='tight',dpi=200)
  ax = data_plot_99.plot(title=par99,marker='+')
  ax.yaxis.grid(True, which='minor') 
  ax.set_xlabel("Number of elements")
  ax.set_ylabel("Parameter value")
  plt.minorticks_on()
  plt.grid(which='minor')
  plt.savefig('N_elements/{}.png'.format(par99),format='png',bbox_inches='tight',dpi=200)
  
  plt.close('all')
  
store.close()