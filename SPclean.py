#####################################################################################################
#
# 
#
# Written by Daniele Michilli
#
#####################################################################################################
import pandas as pd
import tarfile

import time

import RFIexcision
import Group


def openSB(idL,sap,beam):
  #Select the single pulse file
  name = idL+'_SAP'+str(sap)+'_BEAM'+str(beam)
  pulses_file = name+'_singlepulse.tgz'  #per i test

  try:
    pulses_tar = tarfile.open(pulses_file)
    pulses = pulses_tar.extractfile(name+'.singlepulse')

    #Write beam data in a temporary dirty table
    data = pd.read_csv(pulses, delim_whitespace=True)
 
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
    
    data.Sampling = data.Sampling*data.Downfact/2.
    data.rename(columns={'Sampling': 'Down_Time'},inplace=True)

    #Select the interesting columns
    data = data.ix[:,['DM','Sigma','Time','Downfact','Down_Time']]
       
    pulses.close()
    pulses_tar.close()
    
  except IOError:
    a=0  #per i test
    #print "SAP: "+str(sap)+" - BEAM: "+str(beam)+" doesn't exist!"
    data = pd.DataFrame()
    
  return data



def obs_events(idL):
  #Create the table in the memory
  data = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Downfact','Down_Time'])

  #Read each beam data
  for sap in range(0,3):  #0,3
  
    #Incoherent beam
    data_inc = openSB(idL,sap,12)
      
    #########################PER STAMPARE INCOHERENT BEAM
    #data_inc.insert(0,'BEAM',12)
    #data_inc.insert(0,'SAP',sap)
    #data = data.append(data_inc,ignore_index=True)
    ########################   

    for beam in range(13,74):  #13,74
      data_sb = openSB(idL,sap,beam)
      if not data_sb.empty: 
        
        #Remove RFI with sinle-beam techniques
        data_sb = RFIexcision.SB(data_sb)
        
        #Remove RFI with multi-beam techniques
        if not data_inc.empty:
          data_sb = RFIexcision.IB(data_sb,data_inc)
        
        #Insert SAP and BEAM numbers
        data_sb.insert(0,'BEAM',beam)
        data_sb.insert(0,'SAP',sap)
        
        #Add clean data to the table
        data = data.append(data_sb,ignore_index=True)
        
        
  # Si puo diminuire il consumo di memoria facendo partire RFIexcision.MB durante l'elaborazione di sap 3
  #Remove RFI with multi-beam techniques
  data = RFIexcision.MB(data)
  

  
  data = Group.TimeAll(data)  
  
  
  time0=time.clock()
  data = Group.Pulses(data)
  #print 'Time: ',time.clock()-time0,' s'
  
  pulses = Group.Table(data)
  
  #print pulses
  
  #data = RFIexcision.Pulses(data) #,pulses)
  
  
  #Store the table in a DB
  store = pd.HDFStore('SinlgePulses.hdf5','w')
  store[idL] = data
  #store[idL+'_pulses'] = pulses
  store.close()

  return data