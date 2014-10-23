import pandas as pd
import matplotlib as mpl
import numpy as np
import os

#FARE: vedere percentuali assolute e relative con rfi

idL = os.path.basename(os.getcwd())
out_file = open("DB - {}".format(idL),'w')
parameters = open("Parameters - {}".format(idL),'w')

if idL == 'L196709' : DM = 36.35
elif idL == 'L196710' : DM = 43.49
elif idL == 'L202429' : DM = 20.87
elif idL == 'L202433' : DM = 26.18
elif idL == 'L215819' : DM = 43.60
elif idL == 'L94172' : DM = 34.95
elif idL == 'L98710' : DM = 10.68

out_file.write('Observation {}\n\n'.format(idL))
parameters.write('{}\n\n'.format(idL))

if DM < 40.6: step=0.01
elif DM < 141.7: step=0.05
else: step=0.1

DM_hi = DM + step + 0.001
DM_lo = DM - step - 0.001

down_DM_hi = DM_hi - 2.
down_DM_lo = DM_lo - 2.
up_DM_hi = DM_hi + 2.
up_DM_lo = DM_lo + 2.


puls = pd.read_hdf('SinlgePulses.hdf5',idL+'_pulses')
puls = puls[puls.Pulse==0]

rfi_puls = puls[(puls.DM>up_DM_hi) | (puls.DM<down_DM_lo)] 
puls = puls[(puls.DM>DM_lo)&(puls.DM<DM_hi)]

rfi = ( len(puls[(puls.DM>up_DM_lo)&(puls.DM<up_DM_hi)]) + len(puls[(puls.DM>down_DM_lo)&(puls.DM<down_DM_hi)]) ) / 2
num = len(puls)-rfi

data = pd.read_hdf('SinlgePulses.hdf5',idL)

data_all_rfi = data[data.Pulse.isin(rfi_puls.index)]
data_all = data[data.Pulse.isin(puls.index)]
gb_rfi = data_all_rfi.groupby('Pulse',sort=False)
gb = data.groupby('Pulse',sort=False)

N_events = gb.DM.count()
N_events_rfi = gb_rfi.DM.count()

pulsA = puls[N_events<9]
pulsB = puls[(N_events>=9)&(N_events<=12)]
pulsC = puls[N_events>16]

rfiA = rfi_puls[N_events_rfi<9]
rfiB = rfi_puls[(N_events_rfi>=9)&(N_events_rfi<=12)]
rfiC = rfi_puls[N_events_rfi>16]

puls=0
rfi_puls=0


for [puls,rfi_puls] in [[pulsA,rfiA],[pulsB,rfiB],[pulsC,rfiC]]:
  rfi = ( len(puls[(puls.DM>up_DM_lo)&(puls.DM<up_DM_hi)]) + len(puls[(puls.DM>down_DM_lo)&(puls.DM<down_DM_hi)]) ) / 2
  num = len(puls)-rfi
  num_rfi = len(rfi_puls)
  data = data_all[data_all.Pulse.isin(puls.index)]
  gb = data.groupby('Pulse',sort=False)
  Sigma_DM_max = data.Sigma[gb.DM.idxmax()].values
  Sigma_DM_min = data.Sigma[gb.DM.idxmin()].values
  Time_DM_max = data.Time[gb.DM.idxmax()].values
  Time_DM_min = data.Time[gb.DM.idxmin()].values
  DM_min = gb.DM.min()
  Sigma_min = gb.Sigma.min()
  N_events = gb.DM.count()
  
  data_rfi = data_all_rfi[data_all_rfi.Pulse.isin(rfi_puls.index)]
  gb_rfi = data_rfi.groupby('Pulse',sort=False)
  Sigma_DM_max_rfi = data_rfi.Sigma[gb_rfi.DM.idxmax()].values
  Sigma_DM_min_rfi = data_rfi.Sigma[gb_rfi.DM.idxmin()].values
  Time_DM_max_rfi = data_rfi.Time[gb_rfi.DM.idxmax()].values
  Time_DM_min_rfi = data_rfi.Time[gb_rfi.DM.idxmin()].values
  DM_min_rfi = gb_rfi.DM.min()
  Sigma_min_rfi = gb_rfi.Sigma.min()
  N_events_rfi = gb_rfi.DM.count()
  


  out_file.write('Pulses formed by a number of elements between {} e {}\n\n'.format(N_events.min(),N_events.max()))

  out_file.write( 'Astrophysical pulses: {}\nRFIs: {}\n\n'.format(num,len(rfi_puls)))

  if num < np.sqrt(rfi):
    out_file.write( 'Pulsar not found')
    exit()

  selected60 = puls
  selected75 = puls
  selected95 = puls
  selected99 = puls
  selected60_rfi = rfi_puls
  selected75_rfi = rfi_puls
  selected95_rfi = rfi_puls
  selected99_rfi = rfi_puls
  
  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))


    
  out_file.write( "\n\n-------------------------------------\n\npuls.Duration\n")
  parameters.write('puls.Duration\n')
  puls['parameter'] = puls.Duration
  sign = '>'
  
  rfi_puls['parameter'] = rfi_puls.Duration
  
  selected60['parameter'] = selected60.Duration
  selected75['parameter'] = selected75.Duration
  selected95['parameter'] = selected95.Duration
  selected99['parameter'] = selected99.Duration

  selected60_rfi['parameter'] = selected60_rfi.Duration
  selected75_rfi['parameter'] = selected75_rfi.Duration
  selected95_rfi['parameter'] = selected95_rfi.Duration
  selected99_rfi['parameter'] = selected99_rfi.Duration
  
  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100

  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  
  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100

  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")

  


  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")
  
  puls['parameter'] = puls.dDM/(N_events-1)/step
  sign = '>'

  out_file.write( "\n\n-------------------------------------\n\npuls.dDM/(N_events-1)/step\n")  
  parameters.write('puls.dDM/(N_events-1)/step\n')
  rfi_puls['parameter'] = rfi_puls.dDM/(N_events_rfi-1)/step

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  N_events60 = gb60.DM.count()
  N_events75 = gb75.DM.count()
  N_events95 = gb95.DM.count()
  N_events99 = gb99.DM.count()

  selected60['parameter'] = selected60.dDM/(N_events60-1)/step
  selected75['parameter'] = selected75.dDM/(N_events75-1)/step
  selected95['parameter'] = selected95.dDM/(N_events95-1)/step
  selected99['parameter'] = selected99.dDM/(N_events99-1)/step

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  N_events60_rfi = gb60_rfi.DM.count()
  N_events75_rfi = gb75_rfi.DM.count()
  N_events95_rfi = gb95_rfi.DM.count()
  N_events99_rfi = gb99_rfi.DM.count()

  selected60_rfi['parameter'] = selected60_rfi.dDM/(N_events60_rfi-1)/step
  selected75_rfi['parameter'] = selected75_rfi.dDM/(N_events75_rfi-1)/step
  selected95_rfi['parameter'] = selected95_rfi.dDM/(N_events95_rfi-1)/step
  selected99_rfi['parameter'] = selected99_rfi.dDM/(N_events99_rfi-1)/step
  
  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")
  
  
  

  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))

  puls['parameter'] = N_events/puls.Sigma
  sign = '>'

  out_file.write( "\n\n-------------------------------------\n\nN_events/puls.Sigma\n")  
  parameters.write('N_events/puls.Sigma\n')
  rfi_puls['parameter'] = N_events_rfi/rfi_puls.Sigma

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  N_events60 = gb60.DM.count()
  N_events75 = gb75.DM.count()
  N_events95 = gb95.DM.count()
  N_events99 = gb99.DM.count()

  selected60['parameter'] = N_events60/selected60.Sigma
  selected75['parameter'] = N_events75/selected75.Sigma
  selected95['parameter'] = N_events95/selected95.Sigma
  selected99['parameter'] = N_events99/selected99.Sigma


  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  N_events60_rfi = gb60_rfi.DM.count()
  N_events75_rfi = gb75_rfi.DM.count()
  N_events95_rfi = gb95_rfi.DM.count()
  N_events99_rfi = gb99_rfi.DM.count()

  selected60_rfi['parameter'] = N_events60_rfi/selected60_rfi.Sigma
  selected75_rfi['parameter'] = N_events75_rfi/selected75_rfi.Sigma
  selected95_rfi['parameter'] = N_events95_rfi/selected95_rfi.Sigma
  selected99_rfi['parameter'] = N_events99_rfi/selected99_rfi.Sigma

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")




  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")

  puls['parameter'] = abs(puls.DM-puls.DM_c)/puls.dDM
  sign = '>'


  out_file.write( "\n\n-------------------------------------\n\nabs(puls.DM-puls.DM_c)/puls.dDM\n")  
  parameters.write('abs(puls.DM-puls.DM_c)/puls.dDM\n')
  rfi_puls['parameter'] = abs(rfi_puls.DM-rfi_puls.DM_c)/rfi_puls.dDM

  selected60['parameter'] = abs(selected60.DM-selected60.DM_c)/selected60.dDM
  selected75['parameter'] = abs(selected75.DM-selected75.DM_c)/selected75.dDM
  selected95['parameter'] = abs(selected95.DM-selected95.DM_c)/selected95.dDM
  selected99['parameter'] = abs(selected99.DM-selected99.DM_c)/selected99.dDM

  selected60_rfi['parameter'] = abs(selected60_rfi.DM-selected60_rfi.DM_c)/selected60_rfi.dDM
  selected75_rfi['parameter'] = abs(selected75_rfi.DM-selected75_rfi.DM_c)/selected75_rfi.dDM
  selected95_rfi['parameter'] = abs(selected95_rfi.DM-selected95_rfi.DM_c)/selected95_rfi.dDM
  selected99_rfi['parameter'] = abs(selected99_rfi.DM-selected99_rfi.DM_c)/selected99_rfi.dDM

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")
  
  


  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))

  puls['parameter'] = abs(puls.DM-puls.DM_c)/puls.dDM**4
  sign = '>'


  out_file.write( "\n\n-------------------------------------\n\nabs(puls.DM-puls.DM_c)/puls.dDM**4\n")  
  parameters.write('abs(puls.DM-puls.DM_c)/puls.dDM**4\n')
  rfi_puls['parameter'] = abs(rfi_puls.DM-rfi_puls.DM_c)/rfi_puls.dDM**4

  selected60['parameter'] = abs(selected60.DM-selected60.DM_c)/selected60.dDM**4
  selected75['parameter'] = abs(selected75.DM-selected75.DM_c)/selected75.dDM**4
  selected95['parameter'] = abs(selected95.DM-selected95.DM_c)/selected95.dDM**4
  selected99['parameter'] = abs(selected99.DM-selected99.DM_c)/selected99.dDM**4

  selected60_rfi['parameter'] = abs(selected60_rfi.DM-selected60_rfi.DM_c)/selected60_rfi.dDM**4
  selected75_rfi['parameter'] = abs(selected75_rfi.DM-selected75_rfi.DM_c)/selected75_rfi.dDM**4
  selected95_rfi['parameter'] = abs(selected95_rfi.DM-selected95_rfi.DM_c)/selected95_rfi.dDM**4
  selected99_rfi['parameter'] = abs(selected99_rfi.DM-selected99_rfi.DM_c)/selected99_rfi.dDM**4

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")




  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")

  puls['parameter'] = puls.Sigma/Sigma_min
  sign = '<'

  out_file.write( "\n\n-------------------------------------\n\npuls.Sigma/Sigma_min\n")  
  parameters.write('puls.Sigma/Sigma_min\n')
  rfi_puls['parameter'] = rfi_puls.Sigma/Sigma_min_rfi

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  Sigma_min60 = gb60.Sigma.min()
  Sigma_min75 = gb75.Sigma.min()
  Sigma_min95 = gb95.Sigma.min()
  Sigma_min99 = gb99.Sigma.min()

  selected60['parameter'] = selected60.Sigma/Sigma_min60
  selected75['parameter'] = selected75.Sigma/Sigma_min75
  selected95['parameter'] = selected95.Sigma/Sigma_min95
  selected99['parameter'] = selected99.Sigma/Sigma_min99

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  Sigma_min60_rfi = gb60_rfi.Sigma.min()
  Sigma_min75_rfi = gb75_rfi.Sigma.min()
  Sigma_min95_rfi = gb95_rfi.Sigma.min()
  Sigma_min99_rfi = gb99_rfi.Sigma.min()

  selected60_rfi['parameter'] = selected60_rfi.Sigma/Sigma_min60_rfi
  selected75_rfi['parameter'] = selected75_rfi.Sigma/Sigma_min75_rfi
  selected95_rfi['parameter'] = selected95_rfi.Sigma/Sigma_min95_rfi
  selected99_rfi['parameter'] = selected99_rfi.Sigma/Sigma_min99_rfi

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")
  
  
  
  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))

  puls['parameter'] = puls.Sigma/Sigma_min**4
  sign = '<'

  out_file.write( "\n\n-------------------------------------\n\npuls.Sigma/Sigma_min**4\n")
  parameters.write('puls.Sigma/Sigma_min**4\n')
  rfi_puls['parameter'] = rfi_puls.Sigma/Sigma_min_rfi**4

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  Sigma_min60 = gb60.Sigma.min()
  Sigma_min75 = gb75.Sigma.min()
  Sigma_min95 = gb95.Sigma.min()
  Sigma_min99 = gb99.Sigma.min()

  selected60['parameter'] = selected60.Sigma/Sigma_min60**4
  selected75['parameter'] = selected75.Sigma/Sigma_min75**4
  selected95['parameter'] = selected95.Sigma/Sigma_min95**4
  selected99['parameter'] = selected99.Sigma/Sigma_min99**4

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  Sigma_min60_rfi = gb60_rfi.Sigma.min()
  Sigma_min75_rfi = gb75_rfi.Sigma.min()
  Sigma_min95_rfi = gb95_rfi.Sigma.min()
  Sigma_min99_rfi = gb99_rfi.Sigma.min()

  selected60_rfi['parameter'] = selected60_rfi.Sigma/Sigma_min60_rfi**4
  selected75_rfi['parameter'] = selected75_rfi.Sigma/Sigma_min75_rfi**4
  selected95_rfi['parameter'] = selected95_rfi.Sigma/Sigma_min95_rfi**4
  selected99_rfi['parameter'] = selected99_rfi.Sigma/Sigma_min99_rfi**4

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")
  

  

  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")

  puls['parameter'] = abs(Sigma_DM_max-Sigma_DM_min)
  sign = '>'

  out_file.write( "\n\n-------------------------------------\n\nabs(Sigma_DM_max-Sigma_DM_min)\n")  
  parameters.write('abs(Sigma_DM_max-Sigma_DM_min)\n')
  rfi_puls['parameter'] = abs(Sigma_DM_max_rfi-Sigma_DM_min_rfi)

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  Sigma_DM_max60 = data.Sigma[gb60.DM.idxmax()].values
  Sigma_DM_max75 = data.Sigma[gb75.DM.idxmax()].values  
  Sigma_DM_max95 = data.Sigma[gb95.DM.idxmax()].values
  Sigma_DM_max99 = data.Sigma[gb99.DM.idxmax()].values
  Sigma_DM_min60 = data.Sigma[gb60.DM.idxmin()].values
  Sigma_DM_min75 = data.Sigma[gb75.DM.idxmin()].values
  Sigma_DM_min95 = data.Sigma[gb95.DM.idxmin()].values
  Sigma_DM_min99 = data.Sigma[gb99.DM.idxmin()].values

  selected60['parameter'] = abs(Sigma_DM_max60-Sigma_DM_min60)
  selected75['parameter'] = abs(Sigma_DM_max75-Sigma_DM_min75)
  selected95['parameter'] = abs(Sigma_DM_max95-Sigma_DM_min95)
  selected99['parameter'] = abs(Sigma_DM_max99-Sigma_DM_min99)

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  Sigma_DM_max60_rfi = data_rfi.Sigma[gb60_rfi.DM.idxmax()].values
  Sigma_DM_max75_rfi = data_rfi.Sigma[gb75_rfi.DM.idxmax()].values
  Sigma_DM_max95_rfi = data_rfi.Sigma[gb95_rfi.DM.idxmax()].values
  Sigma_DM_max99_rfi = data_rfi.Sigma[gb99_rfi.DM.idxmax()].values
  Sigma_DM_min60_rfi = data_rfi.Sigma[gb60_rfi.DM.idxmin()].values
  Sigma_DM_min75_rfi = data_rfi.Sigma[gb75_rfi.DM.idxmin()].values
  Sigma_DM_min95_rfi = data_rfi.Sigma[gb95_rfi.DM.idxmin()].values
  Sigma_DM_min99_rfi = data_rfi.Sigma[gb99_rfi.DM.idxmin()].values

  selected60_rfi['parameter'] = abs(Sigma_DM_max60_rfi-Sigma_DM_min60_rfi)
  selected75_rfi['parameter'] = abs(Sigma_DM_max75_rfi-Sigma_DM_min75_rfi)
  selected95_rfi['parameter'] = abs(Sigma_DM_max95_rfi-Sigma_DM_min95_rfi)
  selected99_rfi['parameter'] = abs(Sigma_DM_max99_rfi-Sigma_DM_min99_rfi)

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")




  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))

  puls['parameter'] = puls.dDM/puls.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max)
  sign = '<'

  out_file.write( "\n\n-------------------------------------\n\npuls.dDM/puls.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max)\n")  
  parameters.write('puls.dDM/puls.dTime*(Time_DM_min-Time_DM_max)/abs(Time_DM_min-Time_DM_max)\n')
  rfi_puls['parameter'] = rfi_puls.dDM/rfi_puls.dTime*(Time_DM_min_rfi-Time_DM_max_rfi)/abs(Time_DM_min_rfi-Time_DM_max_rfi)

  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)
  Time_DM_max60 = data.Time[gb60.DM.idxmax()].values
  Time_DM_max75 = data.Time[gb75.DM.idxmax()].values  
  Time_DM_max95 = data.Time[gb95.DM.idxmax()].values
  Time_DM_max99 = data.Time[gb99.DM.idxmax()].values
  Time_DM_min60 = data.Time[gb60.DM.idxmin()].values
  Time_DM_min75 = data.Time[gb75.DM.idxmin()].values
  Time_DM_min95 = data.Time[gb95.DM.idxmin()].values
  Time_DM_min99 = data.Time[gb99.DM.idxmin()].values

  selected60['parameter'] = selected60.dDM/selected60.dTime*(Time_DM_min60-Time_DM_max60)/abs(Time_DM_min60-Time_DM_max60)
  selected75['parameter'] = selected75.dDM/selected75.dTime*(Time_DM_min75-Time_DM_max75)/abs(Time_DM_min75-Time_DM_max75)
  selected95['parameter'] = selected95.dDM/selected95.dTime*(Time_DM_min95-Time_DM_max95)/abs(Time_DM_min95-Time_DM_max95)
  selected99['parameter'] = selected99.dDM/selected99.dTime*(Time_DM_min99-Time_DM_max99)/abs(Time_DM_min99-Time_DM_max99)

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  Time_DM_max60_rfi = data_rfi.Time[gb60_rfi.DM.idxmax()].values
  Time_DM_max75_rfi = data_rfi.Time[gb75_rfi.DM.idxmax()].values
  Time_DM_max95_rfi = data_rfi.Time[gb95_rfi.DM.idxmax()].values
  Time_DM_max99_rfi = data_rfi.Time[gb99_rfi.DM.idxmax()].values
  Time_DM_min60_rfi = data_rfi.Time[gb60_rfi.DM.idxmin()].values
  Time_DM_min75_rfi = data_rfi.Time[gb75_rfi.DM.idxmin()].values
  Time_DM_min95_rfi = data_rfi.Time[gb95_rfi.DM.idxmin()].values
  Time_DM_min99_rfi = data_rfi.Time[gb99_rfi.DM.idxmin()].values

  selected60_rfi['parameter'] = selected60_rfi.dDM/selected60_rfi.dTime*(Time_DM_min60_rfi-Time_DM_max60_rfi)/abs(Time_DM_min60_rfi-Time_DM_max60_rfi)
  selected75_rfi['parameter'] = selected75_rfi.dDM/selected75_rfi.dTime*(Time_DM_min75_rfi-Time_DM_max75_rfi)/abs(Time_DM_min75_rfi-Time_DM_max75_rfi)
  selected95_rfi['parameter'] = selected95_rfi.dDM/selected95_rfi.dTime*(Time_DM_min95_rfi-Time_DM_max95_rfi)/abs(Time_DM_min95_rfi-Time_DM_max95_rfi)
  selected99_rfi['parameter'] = selected99_rfi.dDM/selected99_rfi.dTime*(Time_DM_min99_rfi-Time_DM_max99_rfi)/abs(Time_DM_min99_rfi-Time_DM_max99_rfi)


  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")




  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")

  f = data[data.Pulse.isin(puls.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f_rfi = data_rfi[data_rfi.Pulse.isin(rfi_puls.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  if f.empty: f=float('NaN')
  if f_rfi.empty: f_rfi=float('NaN')

  
  puls['parameter'] = f
  sign = '>'

  out_file.write( "\n\n-------------------------------------\n\nf (>)\n")  
  parameters.write('f (>)\n')
  rfi_puls['parameter'] = f_rfi
  
  f60 = data[data.Pulse.isin(selected60.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f75 = data[data.Pulse.isin(selected75.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f95 = data[data.Pulse.isin(selected95.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f99 = data[data.Pulse.isin(selected99.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])

  if f60.empty: f60=float('NaN')
  if f75.empty: f75=float('NaN')
  if f95.empty: f95=float('NaN')
  if f99.empty: f99=float('NaN')  
  
  selected60['parameter'] = f60
  selected75['parameter'] = f75
  selected95['parameter'] = f95
  selected99['parameter'] = f99

  f60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])

  if f60_rfi.empty: f60_rfi=float('NaN')
  if f75_rfi.empty: f75_rfi=float('NaN')
  if f95_rfi.empty: f95_rfi=float('NaN')
  if f99_rfi.empty: f99_rfi=float('NaN')
  
  selected60_rfi['parameter'] = f60_rfi
  selected75_rfi['parameter'] = f75_rfi
  selected95_rfi['parameter'] = f95_rfi
  selected99_rfi['parameter'] = f99_rfi

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")
  



  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))

  puls['parameter'] = f
  sign = '<'

  out_file.write( "\n\n-------------------------------------\n\nf (<)\n")  
  parameters.write('f (<)\n')
  rfi_puls['parameter'] = f_rfi
  
  f60 = data[data.Pulse.isin(selected60.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f75 = data[data.Pulse.isin(selected75.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f95 = data[data.Pulse.isin(selected95.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f99 = data[data.Pulse.isin(selected99.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])

  if f60.empty: f60=float('NaN')
  if f75.empty: f75=float('NaN')
  if f95.empty: f95=float('NaN')
  if f99.empty: f99=float('NaN')  
  
  selected60['parameter'] = f60
  selected75['parameter'] = f75
  selected95['parameter'] = f95
  selected99['parameter'] = f99

  f60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])
  f99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].loc[:,['DM','Time','Pulse']].astype(np.float64).groupby('Pulse',sort=False).apply(lambda x: np.polyfit(x.Time.astype(np.float64),x.DM.astype(np.float64),1)[0])

  if f60_rfi.empty: f60_rfi=float('NaN')
  if f75_rfi.empty: f75_rfi=float('NaN')
  if f95_rfi.empty: f95_rfi=float('NaN')
  if f99_rfi.empty: f99_rfi=float('NaN')
  
  selected60_rfi['parameter'] = f60_rfi
  selected75_rfi['parameter'] = f75_rfi
  selected95_rfi['parameter'] = f95_rfi
  selected99_rfi['parameter'] = f99_rfi

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")




  out_file.write('\nTop candidates:\n60%:{}\n75%:{}\n95%:{}\n99%:{}'.format(selected60[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected75[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected95[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string(),selected99[['Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")

  F_MIN = 119.43  #MHz
  F_MAX = F_MIN + 31.64  #MHz
  k = 4.1488078e3  #s
  delay = np.float32(.6 * k * (F_MIN**-2 - F_MAX**-2))
  data.Time += data.DM * delay
  data_rfi.Time += data_rfi.DM * delay
  puls.Time += puls.DM * delay
  rfi_puls.Time += rfi_puls.DM * delay
  gb = data.groupby('Pulse',sort=False)
  gb_rfi = data_rfi.groupby('Pulse',sort=False)
  
  gb60 = data[data.Pulse.isin(selected60.index)].groupby('Pulse',sort=False)
  gb75 = data[data.Pulse.isin(selected75.index)].groupby('Pulse',sort=False)
  gb95 = data[data.Pulse.isin(selected95.index)].groupby('Pulse',sort=False)
  gb99 = data[data.Pulse.isin(selected99.index)].groupby('Pulse',sort=False)

  gb60_rfi = data_rfi[data_rfi.Pulse.isin(selected60_rfi.index)].groupby('Pulse',sort=False)
  gb75_rfi = data_rfi[data_rfi.Pulse.isin(selected75_rfi.index)].groupby('Pulse',sort=False)
  gb95_rfi = data_rfi[data_rfi.Pulse.isin(selected95_rfi.index)].groupby('Pulse',sort=False)
  gb99_rfi = data_rfi[data_rfi.Pulse.isin(selected99_rfi.index)].groupby('Pulse',sort=False)
  
  out_file.write( "\n\n-------------------------------------\n\ngb.Time.apply(np.var)\n") 
  parameters.write('gb.Time.apply(np.var)\n')
  rfi_puls['parameter'] = gb_rfi.Time.apply(np.var)
  
  puls['parameter'] = gb.Time.apply(np.var)
  sign = '>'

  selected60['parameter'] = gb60.Time.apply(np.var)
  selected75['parameter'] = gb75.Time.apply(np.var)
  selected95['parameter'] = gb95.Time.apply(np.var)
  selected99['parameter'] = gb99.Time.apply(np.var)

  selected60_rfi['parameter'] = gb60_rfi.Time.apply(np.var)
  selected75_rfi['parameter'] = gb75_rfi.Time.apply(np.var)
  selected95_rfi['parameter'] = gb95_rfi.Time.apply(np.var)
  selected99_rfi['parameter'] = gb99_rfi.Time.apply(np.var)

  if sign =='<': 
    puls = puls.sort('parameter',ascending=False)
    rfi_puls = rfi_puls.sort('parameter',ascending=False)
  elif sign =='>': 
    puls = puls.sort('parameter',ascending=True)
    rfi_puls = rfi_puls.sort('parameter',ascending=True)

  p60 = puls.parameter.iloc[ int(.6*num) ]
  p75 = puls.parameter.iloc[ int(.75*num) ]
  p95 = puls.parameter.iloc[ int(.95*num) ]
  p99 = puls.parameter.iloc[ int(.99*num) ]
  
  len60 = selected60.shape[0]
  len75 = selected75.shape[0]
  len95 = selected95.shape[0]
  len99 = selected99.shape[0]

  len60_rfi = selected60_rfi.shape[0]
  len75_rfi = selected75_rfi.shape[0]
  len95_rfi = selected95_rfi.shape[0]
  len99_rfi = selected99_rfi.shape[0]
  
  if sign == '>': 
    selected60 = selected60[selected60.parameter <= p60]
    selected75 = selected75[selected75.parameter <= p75]
    selected95 = selected95[selected95.parameter <= p95]
    selected99 = selected99[selected99.parameter <= p99]
    
    selected60_rfi = selected60_rfi[selected60_rfi.parameter <= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter <= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter <= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter <= p99]
    
  else:
    selected60 = selected60[selected60.parameter >= p60]
    selected75 = selected75[selected75.parameter >= p75]
    selected95 = selected95[selected95.parameter >= p95]
    selected99 = selected99[selected99.parameter >= p99]
        
    selected60_rfi = selected60_rfi[selected60_rfi.parameter >= p60]
    selected75_rfi = selected75_rfi[selected75_rfi.parameter >= p75]
    selected95_rfi = selected95_rfi[selected95_rfi.parameter >= p95]
    selected99_rfi = selected99_rfi[selected99_rfi.parameter >= p99]

  removed60 = ( 1 - (selected60.shape[0]+0.0000001) / float(len60+0.0000001) ) * 100
  removed75 = ( 1 - (selected75.shape[0]+0.0000001) / float(len75+0.0000001) ) * 100
  removed95 = ( 1 - (selected95.shape[0]+0.0000001) / float(len95+0.0000001) ) * 100
  removed99 = ( 1 - (selected99.shape[0]+0.0000001) / float(len99+0.0000001) ) * 100
  
  removed60_rfi = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(len60_rfi+0.0000001) ) * 100
  removed75_rfi = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(len75_rfi+0.0000001) ) * 100
  removed95_rfi = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(len95_rfi+0.0000001) ) * 100
  removed99_rfi = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(len99_rfi+0.0000001) ) * 100
  
  removed60_tot = ( 1 - (selected60.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed75_tot = ( 1 - (selected75.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed95_tot = ( 1 - (selected95.shape[0]+0.0000001) / float(num+0.0000001) ) * 100
  removed99_tot = ( 1 - (selected99.shape[0]+0.0000001) / float(num+0.0000001) ) * 100

  removed60_rfi_tot = ( 1 - (selected60_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed75_rfi_tot = ( 1 - (selected75_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed95_rfi_tot = ( 1 - (selected95_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  removed99_rfi_tot = ( 1 - (selected99_rfi.shape[0]+0.0000001) / float(num_rfi+0.0000001) ) * 100
  
  out_file.write( '60%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p60,removed60,removed60_tot,selected60.shape[0],removed60_rfi,removed60_rfi_tot,selected60_rfi.shape[0]))
  out_file.write( '75%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p75,removed75,removed75_tot,selected75.shape[0],removed75_rfi,removed75_rfi_tot,selected75_rfi.shape[0]))
  out_file.write( '95%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p95,removed95,removed95_tot,selected95.shape[0],removed95_rfi,removed95_rfi_tot,selected95_rfi.shape[0]))
  out_file.write( '99%: {:.2e} - events removed: {:.1f}% ({:.1f}% of tot), {} left - rfi removed: {:.1f}% ({:.1f}% of tot), {} left\n'.format(p99,removed99,removed99_tot,selected99.shape[0],removed99_rfi,removed99_rfi_tot,selected99_rfi.shape[0]))

  parameters.write("{}\n".format(p60))
  parameters.write("{}\n".format(p75))
  parameters.write("{}\n".format(p95))
  parameters.write("{}\n".format(p99))
  parameters.write("\n")


  out_file.write('\nTop candidates:\n{}\n\n\n\n'.format(selected60[['DM','Sigma']].sort('Sigma',ascending=False).head(10).T.to_string()))
  out_file.write("\f")





out_file.write( "\n\n-------------------------------------\n\n\nChannel excision:\n")



count,div = np.histogram(puls.Time,bins=36000)
out_file.write( "36000 bins:\n")
out_file.write( ' 1*median {:.1f}%\n'.format( div.argsort()[count>=1.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 2*median {:.1f}%\n'.format( div.argsort()[count>=2.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 3*median {:.1f}%\n'.format( div.argsort()[count>=3.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))


count,div = np.histogram(puls.Time,bins=3600)
out_file.write( "3600 bins\n")
out_file.write( ' 1*median {:.1f}%\n'.format( div.argsort()[count>=1.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 2*median {:.1f}%\n'.format( div.argsort()[count>=2.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 3*median {:.1f}%\n'.format( div.argsort()[count>=3.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))


count,div = np.histogram(puls.Time,bins=360)
out_file.write( "360 bins\n")
out_file.write( ' 1*median {:.1f}%\n'.format( div.argsort()[count>=1.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 2*median {:.1f}%\n'.format( div.argsort()[count>=2.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))
out_file.write( ' 3*median {:.1f}%\n'.format( div.argsort()[count>=3.*np.median(count[count>0])].shape[0]/float(div.shape[0])*100 ))


out_file.close()