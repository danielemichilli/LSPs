#############################
#
# Single Pulse cleaner
#
# Written by Daniele Michilli
#
#############################

import pandas as pd
import numpy as np
import tarfile

import RFIexcision
import Group
import LSPplot


def openSB(idL,sap,beam):
  #------------------------------------------
  # Creates a table for one .singlepulse file
  #------------------------------------------
  
  name = idL+'_SAP'+str(sap)+'_BEAM'+str(beam)
  pulses_file = name+'_singlepulse.tgz'  #per i test

  try:
    #Open the file
    pulses_tar = tarfile.open(pulses_file)
    pulses = pulses_tar.extractfile(name+'.singlepulse')

    #Create the table
    data = pd.read_csv(pulses, delim_whitespace=True, dtype=np.float32)
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
    data.Sampling = data.Sampling*data.Downfact/2.
    data.rename(columns={'Sampling': 'Duration'},inplace=True)
    data.Duration = data.Duration.astype(np.float32)
    
    #Select the interesting columns
    data = data.ix[:,['DM','Sigma','Time','Duration']]
    
    pulses.close()
    pulses_tar.close()
    
  except IOError:
    #Handle missing beams
    a=0  #per i test
    #print "SAP: "+str(sap)+" - BEAM: "+str(beam)+" doesn't exist!"
    data = pd.DataFrame()
    
  return data



def obs_events(idL):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------
 
  #Create the table in the memory
  data = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse'])
  puls = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse','dDM','DM_min','Sigma_min','dTime','DM_c','Time_c','Sigma_DM_max','Sigma_DM_min','N_events'])

  data = data.astype(np.float32)
  data.SAP = data.SAP.astype(np.uint8)
  data.BEAM = data.BEAM.astype(np.uint8)
  data.Pulse = data.Pulse.astype(np.int32)
  
  puls = puls.astype(np.float32)
  puls.SAP = puls.SAP.astype(np.uint8)
  puls.BEAM = puls.BEAM.astype(np.uint8)
  puls.Pulse = puls.Pulse.astype(np.int8)
  puls.N_events = puls.N_events.astype(np.int16)
  
  
  #Adds each clean beam to the table
  for sap in range(0,3):
  
    #Cleans and groups incoherent beams
    data_inc = openSB(idL,sap,12)
    if not data_inc.empty: 
      data_inc = RFIexcision.IB_Event_Thresh(data_inc) 
      puls_inc = Group.Pulses(data_inc)  
      puls_inc = RFIexcision.IB_Pulse_Thresh(puls_inc)
      puls_inc_copy = puls_inc
        
    #Cleans and groups coherent beams
    for beam in range(13,74):
      data_sb = openSB(idL,sap,beam)
      if not data_sb.empty:
        data_sb = RFIexcision.Event_Thresh(data_sb)
        puls_sb = Group.Pulses(data_sb)
                
        data_sb.insert(0,'BEAM',beam)
        data_sb.insert(0,'SAP',sap)
        data_sb.SAP = data_sb.SAP.astype(np.uint8)
        data_sb.BEAM = data_sb.BEAM.astype(np.uint8)
        puls_sb.insert(0,'BEAM',beam)
        puls_sb.insert(0,'SAP',sap)
        puls_sb.SAP = puls_sb.SAP.astype(np.uint8)
        puls_sb.BEAM = puls_sb.BEAM.astype(np.uint8)
        
        if not data_inc.empty: 
          RFIexcision.Compare_IB(puls_sb,puls_inc,puls_inc_copy)
        
        data = data.append(data_sb,ignore_index=True)        
        puls = puls.append(puls_sb,ignore_index=False)
            
    if not data_inc.empty: 
      data_inc.insert(0,'BEAM',12)
      data_inc.insert(0,'SAP',sap)
      data_inc.SAP = data_inc.SAP.astype(np.uint8)
      data_inc.BEAM = data_inc.BEAM.astype(np.uint8)
      puls_inc.insert(0,'BEAM',12)
      puls_inc.insert(0,'SAP',sap)
      puls_inc.SAP = puls_inc.SAP.astype(np.uint8)
      puls_inc.BEAM = puls_inc.BEAM.astype(np.uint8)
      data = data.append(data_inc,ignore_index=True)
      puls = puls.append(puls_inc,ignore_index=False)
      
  #LSPplot.RFI_channels(data)
  
  #Compares pulses in different beams
  puls = RFIexcision.Compare_Beams(puls[puls.BEAM>12])
    
  #puls = puls.ix[:,['SAP','BEAM','DM','Sigma','Time','Duration','Pulse']]
  
  #data.Time = Group.TimeAlign(data.Time,data.DM)
  #puls.Time = Group.TimeAlign(puls.Time,puls.DM)
  #puls.Time_c = Group.TimeAlign(puls.Time_c,puls.DM)
  

  #Stores the table into a DataBase
  if not data.empty:
    store = pd.HDFStore('SinlgePulses.hdf5','w')
    store.append(idL,data,data_columns=['Pulse'])
    store.append(idL+'_pulses',puls,data_columns=['Pulse'])
    store.close()
    
  return