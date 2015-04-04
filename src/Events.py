import pandas as pd
import numpy as np
import tarfile
import os
import logging

import C_Funct
from Parameters import *

def Loader(folder,idL,sap,beam):
  #------------------------------------------
  # Creates a table for one .singlepulse file
  #------------------------------------------
  
  events = pd.DataFrame()
  meta_data = pd.DataFrame()
      
  name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
  path = 'SAP{}/{}/BEAM{}_sift/sp/'.format(sap,name,beam)
  #events_path = '{}{}/{}{}_singlepulse.tgz'.format(folder,idL,path,name)
  
  events_path = '{}{}/{}_singlepulse.tgz'.format(folder,idL,name)

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
    data['Pulse'] = 0
    data.Pulse = data.Pulse.astype(np.int32)
    
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

  return events,meta_data



def Thresh(events):
  #---------------------------------
  # Applies thresholds to the events
  #---------------------------------

  #Remove low DM events
  cleaned = events[events.DM>DM_MIN]
  
  return cleaned



def Group(events):
  #-----------------------------------
  # Assigns a pulse-code to each event
  #-----------------------------------

  events.sort('DM',inplace=True)
     
  C_Funct.Get_Group(events.DM.values,events.Sigma.values,events.Time.values,events.Duration.values,events.Pulse.values)
  
  events.Pulse = (events.Pulse*np.int32(10)+events.SAP)*np.int32(100)+events.BEAM
  
  return



def TimeAlign(Time,DM):
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 4149. / 2  #s-1
  delay1 = np.float32(k * (F_MIN**-2 - F_MAX**-2))
  
  
  delay2 = np.float32( 0.50833285845243015 * np.digitize(DM, np.arange(2.53*2, 546.48, 2.53)) )
  
  
  #delay2 = delay1 * 2.53 * np.digitize(DM, np.arange(2.53*2, 546.48, 2.53))
  
  #bin=(F_MAX-F_MIN)/288  #288 subbands
  #f_range = F_MAX-F_MIN
  #delay2 = np.float32(((F_MIN+f_range/2)**-2 - (F_MIN+f_range*145/288)**-2) / k)
  
  Time += DM * delay1 + delay2
  
  return Time