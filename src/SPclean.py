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

import Events
import Pulses
import RFIexcision
import LSPplot
import Output
import Candidates
import Internet
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
  dirs = zip(range(saps)*beams,range(beams)*saps)
  
  pulses = pulses_parallel(folder,idL,dirs)
  
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'w')
  for file in os.listdir('{}{}/sp'.format(folder,idL)):
    if file.endswith('.tmp'):
      events = pd.read_hdf('{}{}/sp/{}'.format(folder,idL,file),'events')
      store.append('events',events,data_columns=['Pulse'])
      events = 0
      meta_data = pd.read_hdf('{}{}/sp/{}'.format(folder,idL,file),'meta_data')
      meta_data.reset_index(inplace=True,drop=True)
      store.append('meta_data',meta_data)
      meta_data = 0
      os.remove('{}{}/sp/{}'.format(folder,idL,file))
  store.close()
  
  if pulses.empty: 
    logging.warning("No pulse detected!")
    return
  pulses.sort_index(inplace=True)



  #TO BE TESTED

  '''
  #Remove time-spans affectd by RFI
  pulses.sort_index(inplace=True)
  pulses.Pulse.loc[RFIexcision.time_span(pulses[pulses.BEAM==12])] += 1
  pulses.Pulse.loc[RFIexcision.time_span(pulses[pulses.BEAM>12])] += 1
  pulses = pulses[pulses.Pulse <= RFI_percent]

  #Compare different beams
  pulses.Pulse.loc[RFIexcision.Compare_Beams(pulses[pulses.BEAM>12].copy())] += 1
  pulses = pulses[pulses.Pulse <= RFI_percent]

  #Remove pulses appearing in too many beams
  pulses.Pulse += pulses.apply(lambda x: RFIexcision.puls_beams_select(x,pulses),axis=1).astype(np.int8)
  pulses = pulses[pulses.Pulse <= RFI_percent]
  '''
  
  
  cands = Candidates.candidates(pulses,idL)
  
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  store.append('pulses',pulses,data_columns=['Pulse'])
  if not cands.empty:
    cands.sort(['Rank','Sigma'],ascending=[0,1],inplace=True)
    store.append('candidates',cands,data_columns=['Candidate'])
    cands = cands[cands.main_cand==0].head(30)
  store.close()
    
  if pulses[pulses.Pulse==0].empty: logging.warning("Any reliable pulse detected!")
  else:
    meta_data = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'meta_data')
    
    #Produce the output
    Output.output(folder,idL,pulses,meta_data,cands)
    
  if cands.empty: logging.warning("Any reliable candidate detected!")
  else:
    #Store the best candidates online
    try: Internet.upload(cands,folder,idL)
    except: logging.warning("ATTENTION! Website currently down. Try to upload the observation later with the Upload.py script.")
    Internet.upload_sheet(cands,idL)
    
  return


def pulses_parallel(folder,idL,dirs):
  #Create events, meta_data and pulses lists
  CPUs = mp.cpu_count()
  dirs_range = int(np.ceil(len(dirs)/float(CPUs)))
  
  pool = mp.Pool(CPUs)
  results = [pool.apply_async(lists_creation, args=(folder,idL,dirs[i*dirs_range:(i+1)*dirs_range])) for i in range(CPUs) if len(dirs[i*dirs_range:(i+1)*dirs_range]) > 0]
  pool.close()
  pool.join()
  return pd.concat([p.get() for p in results])


def lists_creation(folder,idL,dirs):
  result = pd.DataFrame()
  for (sap,beam) in dirs:
    
    #Import the events
    events, meta_data = Events.Loader(folder,idL,sap,beam)
    pulses = pd.DataFrame()
    
    if not events.empty:
      try:
        #Correct for the time misalignment of events
        events.sort(['DM','Time'],inplace=True)
        events.Time = Events.TimeAlign(events.Time.copy(),events.DM)
        
        #Group the events
        events.sort(['DM','Time'],inplace=True)
        Events.Group(events)
        
        events.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'events',mode='w')
        meta_data.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'meta_data',mode='a')

        #Apply the thresholds to the events
        events = events[events.Pulse>=0]
        events = Events.Thresh(events)

        #Generate the pulses
        pulses = Pulses.Generator(events)
        pulses = pulses[pulses.Sigma >= 6.5]
        events = events[events.Pulse.isin(pulses.index)]
              
        #Apply global RFI filters to the pulses
        RFIexcision.global_filters(pulses,events)
        pulses = pulses[pulses.Pulse <= RFI_percent]

        #Set a maximum amout of pulses to prevent bad observations to block the pipeline
        pulses.sort('Sigma',ascending=False,inplace=True)
        pulses = pulses.iloc[:3e4]
        events = events[events.Pulse.isin(pulses.index)]
        
        #Apply local RFI filters to the pulses
        RFIexcision.local_filters(pulses,events)
        pulses = pulses[pulses.Pulse <= RFI_percent]
        #events = events[events.Pulse.isin(pulses.index)]
        
        #A set of known pulses is necessary
        #Apply multimoment analysis to the pulses
        #RFIexcision.multimoment(pulses,idL)
        #pulses = pulses[pulses.Pulse <= RFI_percent]
        #events = events[events.Pulse.isin(pulses.index)]
        
      except:
        logging.warning("Some problem arised processing SAP "+str(sap)+" - BEAM "+str(beam)+", it will be discarded")
        pulses = pd.DataFrame()
      
    result = pd.concat((pulses,result))
  
  return result

