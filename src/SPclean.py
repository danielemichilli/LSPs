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
  
  #pulses = lists_creation((folder,idL,2,56))
  
  
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
  
  cands = Candidates.candidates(pulses)
  
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  store.append('pulses',pulses,data_columns=['Pulse'])
  store.append('candidates',cands,data_columns=['Candidate'])
  store.close()
    
  #Compares pulses in different beams
  #pulses.Pulse[pulses.BEAM>12] = 
  RFIexcision.Compare_Beams(pulses[pulses.BEAM>12])
  #pulses.Pulse[pulses.BEAM==12] = 
  RFIexcision.time_span(pulses[pulses.BEAM==12])
  
  #Clean the pulses table
  pulses = pulses[pulses.Pulse <= RFI_percent]
  
  if pulses[pulses.Pulse==0].empty: logging.warning("Any reliable pulse detected!")
  else:
    
    events = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'events',where=['Pulse==pulses.index.tolist()'])
    meta_data = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'meta_data')
    
    #Produce the output
    Output.output(folder,idL,pulses,events,meta_data,cands)
  
  return



def lists_creation((folder,idL,sap,beam)):

  #Import the events
  events, meta_data = Events.Loader(folder,idL,sap,beam)
  pulses = pd.DataFrame()
  
  if not events.empty:
    try:
      #Correct for the time misalignment of events
      events.sort(['DM','Time'],inplace=True)
      events.Time = Events.TimeAlign(events.Time,events.DM)
      
      #Group the events
      events.sort(['DM','Time'],inplace=True)
      Events.Group(events)
      
      events.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'events',mode='w')
      meta_data.to_hdf('{}{}/sp/SAP{}_BEAM{}.tmp'.format(folder,idL,sap,beam),'meta_data',mode='a')

      #Apply the thresholds to the events
      events = Events.Thresh(events)

      #Generate the pulses
      events = events[events.Pulse>=0]
      pulses = Pulses.Generator(events)
      
      #Apply RFI filters to the pulses
      pulses = pulses[pulses.Sigma >= 6.5]
      events = events[events.Pulse.isin(pulses.index)]
      
      RFIexcision.Pulse_Thresh(pulses,events)
      pulses = pulses[pulses.Pulse <= RFI_percent]
      
    except:
      logging.warning("Some problem arised processing SAP "+str(sap)+" - BEAM "+str(beam)+", it will be discarded")
      pulses = pd.DataFrame()
      
  return pulses