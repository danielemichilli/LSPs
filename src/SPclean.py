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
import logging
import math

import Events
import Pulses
import RFIexcision
import LSPplot
from Parameters import *

def obs_events(folder,idL,load_events=False,conf=False):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------

  if conf:
    logging.warning("A confirmation observation will be processed.")
    saps = 1
    beams = 128
  else: 
    saps = 3
    beams = 74
  folders = zip(range(saps)*beams,range(beams)*saps)
  
  #Create events, meta_data and pulses lists
  pool = mp.Pool()
  results = pool.map(lists_creation, [(folder,idL,sap,beam) for (sap,beam) in folders])
  pool.close()
  pool.join()

  pulses = pd.concat(results)
  results = 0
  
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'w')
  for file in os.listdir('{}{}/sp'.format(folder,idL)):
    if file.endswith('.tmp'):
      events = pd.read_hdf('{}{}/sp/{}'.format(folder,idL,file),'events')
      store.append('events',events,data_columns=['Pulse'])
      events = 0
      meta_data = pd.read_hdf('{}{}/sp/{}'.format(folder,idL,file),'meta_data')
      store.append('meta_data',meta_data)
      meta_data = 0
      os.remove('{}{}/sp/{}'.format(folder,idL,file))
  store.close()
  
      
  #store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  #store.append('meta_data',meta_data)
  #store.append('events',events,data_columns=['Pulse'])
  #store.close()
  
  
  if pulses.empty: 
    logging.warning("No pulse detected!")
    return
  pulses.sort_index(inplace=True)

  #cands = candidates(pulses,folders)

  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  store.append('pulses',pulses,data_columns=['Pulse'])
  #store.append('candidates',cands,data_columns=['Candidate'])
  store.close()
    
  #Compares pulses in different beams
  #RFIexcision.Compare_Beams(pulses)

  #Clean the pulses table
  #pulses = pulses[pulses.Pulse <= RFI_percent]
  
  
  
  if pulses[pulses.Pulse==0].empty: logging.warning("Any reliable pulse detected!")
  else:
    
    events = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'events',where=['Pulse==pulses.index.tolist()'])
    meta_data = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'meta_data')
    
    #Produce the output
    pulses = pulses[pulses.Pulse <= 2]
    events = events[events.Pulse.isin(pulses.index)]
    output(folder,idL,pulses,events,meta_data)
  
  return



def lists_creation((folder,idL,sap,beam)):

  #Import the events
  events, meta_data = Events.Loader(folder,idL,sap,beam)
  pulses = pd.DataFrame()
  
  if not events.empty:
    #try:
      #Correct for the time misalignment of events
      events.sort(['DM','Time'],inplace=True)
      events.Time = Events.TimeAlign(events.Time,events.DM)
      
      #Group the events
      events.sort(['DM','Time'],inplace=True)
      Events.Group(events)
      
      events.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'events',mode='w')
      meta_data.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'meta_data',mode='a')
      
      #store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
      #store.append('meta_data',meta_data)
      #store.append('events',events,data_columns=['Pulse'])
      #store.close()

      #Apply the thresholds to the events
      events = Events.Thresh(events)

      #Generate the pulses
      events = events[events.Pulse>=0]
      pulses = Pulses.Generator(events)

      #Apply RFI filters to the pulses
      pulses = pulses[pulses.Sigma >= 6.5]
      events = events[events.Pulse.isin(pulses.index)]
      RFIexcision.Pulse_Thresh(pulses,events)

      events = 0
      pulses = pulses[pulses.Pulse < RFI_percent]

    #except:
      logging.warning("Some problem arised processing SAP "+str(sap)+" - BEAM "+str(beam)+", it will be discarded")

  return pulses



  


def output(folder,idL,pulses,events,meta_data):
  pulses.sort('Sigma',ascending=False,inplace=True)
  
  store = '{}{}/sp/beam_plots'.format(folder,idL)
  
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  
  for n in gb_puls.indices.iterkeys():
    name = 'SAP{}_BEAM{}'.format(n[0],n[1])
    os.makedirs('{}/{}'.format(store,name))
    
  #pool = mp.Pool()
  #pool.map(output_beams, [(pulses,events,meta_data,store,idL,n) for n in gb_puls.indices.iterkeys()])
  #pool.close()
  #pool.join()  
  
  for n in gb_puls.indices.iterkeys():
    output_beams((pulses,events,meta_data,store,idL,n))
  
  pulses = pulses[pulses.Pulse==0]
  output_pointing(pulses,folder,idL)
  
  return



def output_beams((pulses,events,meta_data,folder,obs,(sap,beam))):
  
  top = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==0)]
  good = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==1)]
  
  if beam == 12:
    LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates(0-9).png'.format(folder,sap,beam),obs)
    LSPplot.sp_shape(top.iloc[10:20],events,'{}/SAP{}_BEAM{}/top_candidates(10-19).png'.format(folder,sap,beam),obs)
    LSPplot.sp_shape(top.iloc[20:30],events,'{}/SAP{}_BEAM{}/top_candidates(20-29).png'.format(folder,sap,beam),obs)
    LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
    
  else:
    LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates.png'.format(folder,sap,beam),obs)
    LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
  
  return
  
  

def output_pointing(pulses,folder,idL):
  
  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index = b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}{}/sp/files/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',\
    columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],\
    header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  puls = pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30)
  if not puls.empty: LSPplot.obs_top_candidates(puls,store='{}{}/sp/files/inc_top_candidates.png'.format(folder,idL),incoherent=True)
  
  for sap in pulses.SAP.unique():
    puls = pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10)
    if not puls.empty: LSPplot.obs_top_candidates(puls,store='{}{}/sp/files/top_candidates(SAP{}).png'.format(folder,idL,sap))

  return






def candidates(pulses,folders):
  #Create candidates lists
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  pool = mp.Pool()
  results = pool.map(Pulses.Repeated_candidates, [(pulses,sap,beam) for (sap,beam) in folders])
  pool.close()
  pool.join()
  #Pulses.Repeated_candidates(pulses)
  
  pulses.Candidate = pd.concat(results)
  results = 0
  
  #Unify the same repeated candidates in different beams
  pulses.sort('DM',inplace=True)
  for sap in range(3):
    cands_SAP = pulses[(pulses.SAP==sap)&(pulses.Candidate>0)]
    diff_DM = np.abs(cands_SAP.DM-cands_SAP.DM.shift())
    diff_DM.iloc[0] = diff_DM.iloc[1]
    diff_DM -= 0.5
    diff_DM = diff_DM.round()
    cands_SAP = cands_SAP.Candidate
    cands_SAP[diff_DM<.1] = np.nan
    cands_SAP.fillna(method='pad',inplace=True)
    pulses.Candidate.loc[cands_SAP.index] = cands_SAP
    
    
  
  
  
  #si puo' aggiungere una funzione per non avere lo stesso unique candidate in beams differenti
  pulses.sort('Sigma',ascending=False,inplace=True)
  cands_unique = pulses[(pulses.Pulse==0)&(pulses.Candidate==-1)&(pulses.Sigma>=10)].groupby('N_events')['Sigma'].nlargest(4)
  pulses.Candidate.loc[cands_unique.index.get_level_values(1)] = np.arange(cands_unique.shape[0])
  
  cands = candidates_generator(pulses[pulses.Candidate>=0])
  
  return cands


def candidates_generator(pulses):
  #Prende max N_pulses e mean DM tra tutti i beams, non so se sia l'opzione migliore
  cands = pulses.groupby(['Candidate','BEAM']).agg({'Sigma':np.sum,'SAP':np.size,'DM':np.mean}).groupby(level=0).agg({'Sigma':np.max,'SAP':np.max,'DM':np.mean})
  cands.index.name = None
  cands.rename(columns={'SAP': 'N_pulses'}, inplace=True)


