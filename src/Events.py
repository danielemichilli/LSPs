import pandas as pd
import numpy as np
import tarfile
import os
import logging

import C_Funct
from Parameters import *

def Loader(directory,sap,beam):
  #------------------------------------------
  # Creates a table for one .singlepulse file
  #------------------------------------------
  
  events = pd.DataFrame()
  meta_data = pd.DataFrame()
      
  name = os.path.basename(directory).split('_single')[0]
  
  try:
    #Open the file
    tar_file = tarfile.open(directory)
    events_file = tar_file.extractfile(name+'.singlepulse')
    inf_file = tar_file.extractfile(name+'.inf')
    
    #Load the events and meta-data tables
    data = pd.read_csv(events_file, delim_whitespace=True, dtype=np.float32)
    inf = pd.read_csv(inf_file, sep="=", dtype=str,error_bad_lines=False,warn_bad_lines=False,header=None,skipinitialspace=True)

    tar_file.close()
    events_file.close()
    
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
        
    data['Duration'] = data.Sampling*data.Downfact
    data.Duration = data.Duration.astype(np.float32)
    data.Sample = data.Sample.astype(np.int32)
    
    data = data.ix[:,['DM','Sigma','Time','Duration','Sample','Downfact']]
    data.index.name = 'idx'
    
    data.insert(0,'BEAM',beam)
    data.insert(0,'SAP',sap)
    data.SAP = data.SAP.astype(np.uint8)
    data.BEAM = data.BEAM.astype(np.uint8)      
    data['Pulse'] = 0
    data.Pulse = data.Pulse.astype(np.int64)
    
    inf = inf.iloc[[0,1,2,4,5,7],1]
    inf.iloc[0] = inf.iloc[0].replace("_rfifind","")
    inf = pd.DataFrame(inf).T
    inf.columns=['File','Telescope','Instrument','RA','DEC','Epoch']
    inf = inf.astype(str)
    inf['File'] = inf.File.apply(lambda x: os.path.basename(x))
    
    inf.insert(0,'BEAM',beam)
    inf.insert(0,'SAP',sap)
    inf.SAP = inf.SAP.astype(np.uint8)
    inf.BEAM = inf.BEAM.astype(np.uint8)
    
    #Append to the existing tables
    events = data
    meta_data = inf
    
  except (IOError,pd.parser.CParserError):
    #Handle missing beams
    logging.warning("SAP "+str(sap)+" - BEAM "+str(beam)+" doesn't exist")

  return events,meta_data



def Thresh(events):
  #---------------------------------
  # Applies thresholds to the events
  #---------------------------------

  #Remove events at the end of the observation
  events = events[events.Time < DURATION - 10]
  
  #Remove low-DM events
  #events = events[events.DM > DM_MIN]
  
  return events



def Group(events):
  #-----------------------------------
  # Assigns a pulse-code to each event
  #-----------------------------------

  C_Funct.Get_Group(events.DM.values, events.Sigma.values, events.Time.values, events.Duration.values, events.Pulse.values)
  
  events.Pulse = (events.Pulse * 10 + events.SAP) * 1000 + events.BEAM
  
  return



def TimeAlign(Time,DM):
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 4148.808 #s-1
  delay = k * (F_MIN**-2 - F_MAX**-2)
  Time += np.float32( delay * DM / 2 )
  
  return Time

