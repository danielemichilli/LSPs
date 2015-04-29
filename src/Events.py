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

  C_Funct.Get_Group(events.DM.values,events.Sigma.values,events.Time.values,events.Duration.values,events.Pulse.values)
  
  events.Pulse = (events.Pulse*np.int32(10)+events.SAP)*np.int32(100)+events.BEAM
  
  return



def TimeAlign(Time,DM):
  #-------------------------------------------------
  # Corrects for the time misalignment of the pulses
  #-------------------------------------------------
  
  # Quantifies the misalignment for a broad-band pulse
  # Only the extreme frequencies are taken into account
  k = 4149. #s-1
  delay = k * (F_MIN**-2 - F_MAX**-2)

  DM_steps = np.array((2.52,5.05,7.58,10.11,12.64,15.17,17.7,20.23,22.76,25.29,27.82,30.35,32.88,35.41,37.94,40.47,65.81,91.11,116.41,141.71,242.96,344.16,445.36,546.56))

  DM_n = np.digitize(DM, DM_steps) - 1
    
  a = 0.253 * DM_n[DM_n<15]
  a[a==0.253*9] *= 2
  b = 0.253 * 14 + 25.3 * delay * ( DM_n[(DM_n>=15)&(DM_n<19)] - 14 )
  c = 0.253 * 14 + 25.3 * delay * 4 + 4 * 25.3 * delay * ( DM_n[DM_n>=19] - 18 )
  
  DM_n = np.concatenate((a,b,c))
  
  Time += np.float32( delay * DM / 2 + DM_n )
  
  return Time

