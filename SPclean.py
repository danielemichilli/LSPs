#############################
#
# Single Pulse cleaner
#
# Written by Daniele Michilli
#
#############################

import pandas as pd
import tarfile

import RFIexcision
import Group


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
    data = pd.read_csv(pulses, delim_whitespace=True)
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
    data.Sampling = data.Sampling*data.Downfact/2.
    data.rename(columns={'Sampling': 'Duration'},inplace=True)

    #Select the interesting columns
    data = data.ix[:,['DM','Sigma','Time','Downfact','Duration']]
       
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
  data = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Downfact','Duration'])
  puls = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Downfact','Duration'])

  #Adds each clean beam to the table
  for sap in range(0,3):
  
    #Cleans and groups incoherent beams
    data_inc = openSB(idL,sap,12)
    if not data_inc.empty: 
      data_inc = RFIexcision.IB_Event_Thresh(data_inc,puls_inc)
      data_inc = Group.TimeAlign(data_inc) 
      data_inc = Group.Pulses(data_inc)
      puls_inc = Group.Table(data_inc)      
      data_inc = RFIexcision.IB_Pulse_Thresh(data_inc,puls_inc)
      data_inc.insert(0,'BEAM',12)
      data_inc.insert(0,'SAP',sap)
      data = data.append(data_inc,ignore_index=True)
      puls = puls.append(puls_inc,ignore_index=False)
    
    #Cleans and groups coherent beams
    for beam in range(13,74):
      data_sb = openSB(idL,sap,beam)
      if not data_sb.empty: 
        data_sb = RFIexcision.Event_Thresh(data_sb)
        data_sb = Group.TimeAlign(data_sb) 
        data_sb = Group.Pulses(data_sb)
        puls_sb = Group.Table(data_sb)
#        data_sb = RFIexcision.Pulse_Thresh(data_sb,puls_sb)
#        data_sb = RFIexcision.Compare_IB(data_sb,puls_sb)
        data_sb.insert(0,'BEAM',beam)
        data_sb.insert(0,'SAP',sap)
        data = data.append(data_sb,ignore_index=True)        
        puls = puls.append(puls_sb,ignore_index=False)
  
  #Compares pulses in different beams
#  data = RFIexcision.Compare_Beams(data,puls)
  
  
  #Stores the table into a DataBase
  store = pd.HDFStore('SinlgePulses.hdf5','w')
  store[idL] = data
  store[idL+'_pulses'] = puls
  store.close()

  return