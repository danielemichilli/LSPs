#############################
#
# Single Pulse cleaner
#
# Written by Daniele Michilli
#
#############################

import pandas as pd
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')
import multiprocessing as mp
import matplotlib.pyplot as plt

import Events
import Pulses
import RFIexcision
import LSPplot
from Parameters import *

def obs_events(folder,idL):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------

  folders = zip(range(3)*74,range(74)*3)
  
  #folders = [(0,0),(0,68)]
  
  #Create events, meta_data and pulses lists
  #pool = mp.Pool(mp.cpu_count()-1)
  results = [lists_creation((folder,idL,sap,beam)) for (sap,beam) in folders]
  #pool.close()
  #pool.join()
  
  meta_data = pd.DataFrame()
  events = pd.DataFrame()
  pulses = pd.DataFrame()
  
  for idx, (meta,event,puls) in enumerate(results):
    meta_data = meta_data.append(meta,ignore_index=True)
    events = events.append(event,ignore_index=True)
    pulses = pulses.append(puls)
    results[idx] = 0
  results = 0

  #Compares pulses in different beams
  #RFIexcision.Compare_Beams(pulses)

  #Clean the pulses table
  #pulses = pulses[pulses.Pulse <= RFI_percent]

  #Clean the events table
  #events = events.loc[events.Pulse.isin(pulses.index)]

  #Store the events
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'w')
  store.append(idL,events,data_columns=['Pulse'])
  store.close() 



  
  #def ALERT(sap,beam,dm,limit,counts,s_max,s_mean):
    #file = open('{}{}/sp/ALERTS'.format(folder,idL),'a+')
    #if not file.readlines():
      #file.write('Pulsar candidates\n\n')
      #file.write('SAP\tBEAM\tDM\tLimit\tCounts\tS_Max\tS_Mean\n')
    #file.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4}\t{5:.2f}\t{6:.2f}\n'.format(sap,beam,dm,limit,counts,s_max,s_mean))
    #file.close()
    #return
  
  
  
  
  #data = pulses.loc[pulses.BEAM>12,['SAP','BEAM','DM','Sigma','Pulse']]
  #data.DM = data.DM.round(decimals=1)
  #noise_level = data.groupby(['SAP','BEAM','DM'],sort=False).Pulse.count()
  #noise_level = noise_level.mean(level=['SAP','DM']) + 5.* noise_level.std(level=['SAP','DM'])   #5 sigmas tollerance
  #noise_level.dropna(inplace=True)
    
  #beams = data.groupby(['SAP','BEAM'],sort=False)
  #for ind,beam in beams:
    #counts = beam.groupby(['SAP','DM'],sort=False).Pulse.count()
    #counts = counts[(counts>noise_level.loc[counts.index])&(counts>5)]
    #if not counts.empty:
      #for i in counts.index.get_level_values('DM'):
        #Sigma = beams.get_group((ind[0],ind[1]))
        #Sigma = Sigma.Sigma[Sigma.DM==i]
        #ALERT(ind[0],ind[1],i,noise_level.loc[ind[0],i],counts.loc[ind[0],i],Sigma.max(),Sigma.mean())
      



        
        
  
  # ALERT su incoherent beams rispetto a cosa? media dei tre? media di piu' osservazioni?
  
  #noise_level = pulses.loc[pulses.BEAM>12].groupby(['SAP','BEAM','DM'],sort=False).count().Pulse
  #noise_level = noise_level.mean(level='DM')# + 3.* noise_level.std(level='DM')   #3 sigmas tollerance
  #noise_level.dropna(inplace=True)
  
  #beams = pulses.loc[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False)
  #for ind,beam in beams:
    #counts = beam.groupby('DM',sort=False).Pulse.count()
    #counts = counts[counts>noise_level.loc[counts.index]]
    #if not counts.empty:
      #ALERT(beam.SAP.iloc[0],beam.BEAM.iloc[0],counts.index.values)


  #Store the pulses
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  store.append(idL+'_pulses',pulses,data_columns=['Pulse'])
  store.append('meta_data',meta_data)
  store.close()
  
  #Produce the output
  output(folder,idL,pulses,events,meta_data)
  
  return




def lists_creation((folder,idL,sap,beam)):
  #folder = parameters[0]
  #idL = parameters[1]
  #sap = parameters[2]
  #beam = parameters[3]

  #Import the events
  events, meta_data = Events.Loader(folder,idL,sap,beam)
  pulses = pd.DataFrame()

  if not events.empty:
    
    #Apply the thresholds to the events
    events = Events.Thresh(events)
    
    #Correct for the time misalignment of events
    events.sort('DM',inplace=True)
    events.Time = Events.TimeAlign(events.Time,events.DM)
    
    #Group the events
    Events.Group(events)
    #events = events[events.Pulse>0]
    
    #Generate the pulses
    pulses = Pulses.Generator(events)
    
    #Clean the events table
    #events = events.loc[events.Pulse.isin(pulses.index)]

  return (meta_data,events,pulses)




def output(folder,idL,pulses,events,meta_data):
  
  store = '{}{}'.format(folder,idL)
  
  pulses.sort('Sigma',ascending=False,inplace=True)
  
  gb_rfi = pulses.loc[pulses.Pulse==1].groupby(['SAP','BEAM'],sort=False)
  
  pulses = pulses.loc[pulses.Pulse==0]

  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  
  for n in gb_puls.indices.iterkeys():
    name = 'SAP{}_BEAM{}'.format(n[0],n[1])
    os.makedirs('{}/sp/files'.format(store)+name)

  gb_event = events.loc[events.Pulse.isin(gb_puls.head(30).index)].groupby(['SAP','BEAM'],sort=False)
  events = 0  
  
  gb_md = meta_data.groupby(['SAP','BEAM'],sort=False)
  meta_data = 0
  
      
  def Group_Clean(gb,n):
    try:
      return gb.get_group(n)
    except:
      return pd.DataFrame()
  
  
  pool = mp.Pool(mp.cpu_count()-1)
  pool.map(LSPplot.plot, [(gb_puls.get_group(n),gb_rfi.get_group(n),gb_md.get_group(n),gb_event.get_group(n),store) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()
  
  
  #for n in gb_puls.indices.iterkeys():
    #LSPplot.plot((gb_puls.get_group(n),gb_rfi.get_group(n),gb_md.get_group(n),gb_event.get_group(n),store))
  

  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index=b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}{}/sp/files/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',\
    columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],\
    header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  
  LSPplot.obs_top_candidates(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),\
                             store='{}{}/sp/files/inc_top_candidates.png'.format(folder,idL),incoherent=True)

  for sap in pulses.SAP.unique():
    LSPplot.obs_top_candidates(pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10),\
                               store='{}{}/sp/files/top_candidates(SAP{}).png'.format(folder,idL,sap))
    

  pulses = 0

  return
