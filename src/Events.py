import pandas as pd
import numpy as np
import tarfile
import os
import logging

import C_Funct
from Parameters import *

def Loader(parameters):
  #------------------------------------------
  # Creates a table for one .singlepulse file
  #------------------------------------------
  
  folder = parameters[0]
  idL = parameters[1]
  beam = parameters[2]
  
  events = pd.DataFrame()
  meta_data = pd.DataFrame()
  
  for sap in range(0,3):
      
    name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
    path = 'SAP{}/{}/BEAM{}_sift/sp/'.format(sap,name,beam) #'' per i test
    events_path = '{}{}/{}{}_singlepulse.tgz'.format(folder,idL,path,name)
    
    try:
      #Open the file
      tar_file = tarfile.open(events_path)
      events_file = tar_file.extractfile(name+'.singlepulse')
      inf_file = tar_file.extractfile(name+'.inf')
      
      #Load the events and meta-data tables
      data = pd.read_csv(events_file, delim_whitespace=True, dtype=np.float32)
      inf = pd.read_csv(inf_file, sep="=", dtype=str,error_bad_lines=False,warn_bad_lines=False,header=None,skipinitialspace=True)
      
      tar_file.close()
      events_file.close()
      
      data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
      
      data.Sampling = data.Sampling*data.Downfact
      data.rename(columns={'Sampling': 'Duration'},inplace=True)
      data.Duration = data.Duration.astype(np.float32)
      
      data = data.ix[:,['DM','Sigma','Time','Duration']]
      
      data.insert(0,'BEAM',beam)
      data.insert(0,'SAP',sap)
      data.SAP = data.SAP.astype(np.uint8)
      data.BEAM = data.BEAM.astype(np.uint8)      
      
      inf = inf.iloc[[0,1,2,4,5,7],1]
      inf.iloc[0] = inf.iloc[0].replace("_rfifind","")
      inf = pd.DataFrame(inf).T
      inf.columns=['File','Telescope','Instrument','RA','DEC','Epoch']
      inf = inf.astype(str)
      
      inf.insert(0,'BEAM',beam)
      inf.insert(0,'SAP',sap)
      inf.SAP = inf.SAP.astype(np.uint8)
      inf.BEAM = inf.BEAM.astype(np.uint8)
      
      #Append to the existing tables
      events = events.append(data,ignore_index=True)
      meta_data = meta_data.append(inf,ignore_index=False)
      
    except (IOError,pd.parser.CParserError):
      #Handle missing beams
      logging.warning("SAP "+str(sap)+" - BEAM "+str(beam)+" doesn't exist!")

  return (events,meta_data)



def Thresh(events):
  #-----------------------------------------------------
  # Applies thresholds to the events in a coherent beams
  #-----------------------------------------------------

  #Remove low DM events
  events = events[events.DM>DM_MIN]
  
  return events



def Group(events):
  #-----------------------------------
  # Assigns a pulse-code to each event
  #-----------------------------------

  events.sort(['SAP','BEAM','DM'],inplace=True)  
  
  C_Funct.Get_Group(events.DM.values,events.Sigma.values,events.Time.values,events.Duration.values,events.Pulse.values)

  events.Pulse = (events.Pulse*np.int32(10)+events.SAP)*np.int32(100)+events.BEAM
  
  return