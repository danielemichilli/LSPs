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
    inf_file = pulses_tar.extractfile(name+'.inf')
    
    #Create the table
    data = pd.read_csv(pulses, delim_whitespace=True, dtype=np.float32)
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
    data.Sampling = data.Sampling*data.Downfact/2.
    data.rename(columns={'Sampling': 'Duration'},inplace=True)
    data.Duration = data.Duration.astype(np.float32)
    inf = pd.read_csv(inf_file, sep="=", dtype=str,error_bad_lines=False,header=None,skipinitialspace=True)
    
    #Select the interesting columns
    data = data.ix[:,['DM','Sigma','Time','Duration']]
    inf = inf.iloc[[2,3,0,1,4,5,7],1]  #CAMBIARE! eliminare 2 e 3
    inf.iloc[0] = inf.iloc[0].replace("_rfifind","")
    
    pulses.close()
    pulses_tar.close()
    
  except IOError:
    #Handle missing beams
    data = pd.DataFrame()
    inf = pd.DataFrame()
    
  return data,inf



def obs_events(idL,logging):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------
 
  #Creates the tables
  meta_data = pd.DataFrame()
  data = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse'])
  puls = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse','dDM','dTime','DM_c','Time_c'])
  
  data = data.astype(np.float32)
  data.SAP = data.SAP.astype(np.uint8)
  data.BEAM = data.BEAM.astype(np.uint8)
  data.Pulse = data.Pulse.astype(np.int32)
  
  puls = puls.astype(np.float32)
  puls.SAP = puls.SAP.astype(np.uint8)
  puls.BEAM = puls.BEAM.astype(np.uint8)
  puls.Pulse = puls.Pulse.astype(np.int8)
  
  
  #Adds each clean beam to the table
  for sap in range(0,3):
    
    print 'SAP: ',sap
  
    #Cleans and groups incoherent beams
    data_inc,inf = openSB(idL,sap,12)
    if data_inc.empty: logging.info("SAP "+str(sap)+" - BEAM 12 doesn't exist!")
    else:
      data_inc = RFIexcision.IB_Event_Thresh(data_inc) 
      data_inc, puls_inc = Group.Pulses(data_inc,sap,12)
     
      inf.put(3,beam)  #CAMBIARE! inserire i valori invece di put
      inf.put(2,sap)
      meta_data = meta_data.append(inf,ignore_index=False)
        
    #Cleans and groups coherent beams
    for beam in range(13,74):
      data_sb,inf = openSB(idL,sap,beam)
      if data_sb.empty: logging.info("SAP "+str(sap)+" - BEAM "+str(beam)+" doesn't exist!")
      else:
        data_sb = RFIexcision.Event_Thresh(data_sb)
        data_sb, puls_sb = Group.Pulses(data_sb,sap,beam)
        if not data_inc.empty: 
          RFIexcision.Compare_IB(puls_sb,puls_inc)
        
        data_sb.insert(0,'BEAM',beam)
        data_sb.insert(0,'SAP',sap)
        data_sb.SAP = data_sb.SAP.astype(np.uint8)
        data_sb.BEAM = data_sb.BEAM.astype(np.uint8)
        puls_sb.insert(0,'BEAM',beam)
        puls_sb.insert(0,'SAP',sap)
        puls_sb.SAP = puls_sb.SAP.astype(np.uint8)
        puls_sb.BEAM = puls_sb.BEAM.astype(np.uint8)
        data = data.append(data_sb,ignore_index=True)        
        puls = puls.append(puls_sb,ignore_index=False)
        inf.put(3,beam)  #cambiare dtype in int
        inf.put(2,sap)
        print inf
        meta_data = meta_data.append(inf,ignore_index=False)
            
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
      
  #Compares pulses in different beams
  puls = RFIexcision.Compare_Beams(puls[puls.BEAM>12])
  meta_data.columns=['SAP','BEAM','File','Telescope','RA','DEC','Epoch']
  meta_data = meta_data.astype(str)

  #output(puls,data,meta_data)

  #Stores the table into a DataBase
  if not data.empty:
    store = pd.HDFStore('SinlgePulses.hdf5','w')
    store.append(idL,data,data_columns=['Pulse'])
    store.append(idL+'_pulses',puls,data_columns=['Pulse'])
    store.append('meta_data',meta_data)
    store.close()
      
  return


def output(puls,data,meta_data):
  
  puls = puls.sort_by(['SAP','BEAM','Pulse','Sigma'])
  
  for sap in range(0,3):
    for beam in range(12,74):
      puls_plot = puls[(puls.SAP==sap)&(puls.BEAM==beam)]
      data_plot = data[data.Pulse.isin(puls_plot.index)]
      meta_data_plot = meta_data_plot[(meta_data.SAP==sap)&(meta_data.BEAM==beam)]
      if not puls.empty:
        Pulse_min = puls_plot.Pulse.unique().min()
        top_candidates = puls_plot.iloc[:10]
        LSPplot(idL,puls_plot[puls_plot.Pulse==Pulse_min].iloc[10:],puls_plot[puls_plot.Pulse==Pulse_min+1],top_candidates,meta_data_plot,color=True,size=True,store=True)
        
    
        
  return
