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
import os
import logging
import matplotlib as mpl
mpl.use('Agg')

import RFIexcision
import Group
import LSPplot


def initialize():
  #Creates the tables in memory
  meta_data = pd.DataFrame(columns=['SAP','BEAM','File','Telescope','Instrument','RA','DEC','Epoch'])
  data = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse'])
  puls = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse','dDM','dTime','DM_c','Time_c','N_events'])
  
  data = data.astype(np.float32)
  data.SAP = data.SAP.astype(np.uint8)
  data.BEAM = data.BEAM.astype(np.uint8)
  data.Pulse = data.Pulse.astype(np.int32)
  
  puls = puls.astype(np.float32)
  puls.SAP = puls.SAP.astype(np.uint8)
  puls.BEAM = puls.BEAM.astype(np.uint8)
  puls.Pulse = puls.Pulse.astype(np.int8)
  puls.N_events = puls.N_events.astype(np.int16)
  
  return data,puls,meta_data


def openSB(folder,idL,sap,beam):
  #------------------------------------------
  # Creates a table for one .singlepulse file
  #------------------------------------------
  
  name = '{}_SAP{}_BEAM{}'.format(idL,sap,beam)
  path = 'SAP{}/{}/BEAM{}_sift/sp/'.format(sap,name,beam) #'' per i test
  pulses_file = '{}{}/{}{}_singlepulse.tgz'.format(folder,idL,path,name)
  
  try:
    #Open the file
    pulses_tar = tarfile.open(pulses_file)
    pulses = pulses_tar.extractfile(name+'.singlepulse')
    inf_file = pulses_tar.extractfile(name+'.inf')
    
    #Create the table
    data = pd.read_csv(pulses, delim_whitespace=True, dtype=np.float32)
    data.columns = ['DM','Sigma','Time','Sample','Downfact','Sampling','a','b','c']
        
    data.Sampling = data.Sampling*data.Downfact
    data.rename(columns={'Sampling': 'Duration'},inplace=True)
    data.Duration = data.Duration.astype(np.float32)
    inf = pd.read_csv(inf_file, sep="=", dtype=str,error_bad_lines=False,warn_bad_lines=False,header=None,skipinitialspace=True)
    
    #Select the interesting columns
    data = data.ix[:,['DM','Sigma','Time','Duration']]
    inf = inf.iloc[[0,1,2,4,5,7],1]
    inf.iloc[0] = inf.iloc[0].replace("_rfifind","")
    inf = pd.DataFrame(inf).T
    inf.columns=['File','Telescope','Instrument','RA','DEC','Epoch']

    pulses.close()
    pulses_tar.close()
    
  except (IOError,pd.parser.CParserError):
    #Handle missing beams
    data = pd.DataFrame()
    inf = pd.DataFrame()
    
  return data,inf



def obs_events(folder,idL):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------
  
  data,puls,meta_data = initialize()
  puls_inc_all = pd.DataFrame()
  puls_inc = pd.DataFrame()
  
  #Adds each clean beam to the table
  for sap in range(0,3):
    logging.info('SAP: %s',sap)
    
    #Cleans and groups incoherent beams
    data_inc,inf = openSB(folder,idL,sap,12)
    if data_inc.empty: logging.warning("SAP "+str(sap)+" - BEAM 12 doesn't exist!")
    else:
      data_inc = RFIexcision.Event_Thresh(data_inc)
      data_inc, puls_inc = Group.Pulses(data_inc,sap,12)
      
      inf.insert(0,'BEAM','12')
      inf.insert(0,'SAP',str(sap))
      meta_data = meta_data.append(inf,ignore_index=False)
        
    #Cleans and groups coherent beams
    for beam in range(13,74):  #13,74
      data_sb,inf = openSB(folder,idL,sap,beam)
      if data_sb.empty: logging.warning("SAP "+str(sap)+" - BEAM "+str(beam)+" doesn't exist!")
      else:
        data_sb = RFIexcision.Event_Thresh(data_sb)
        data_sb, puls_sb = Group.Pulses(data_sb,sap,beam)
        #if not data_inc.empty: 
          #puls_sb,puls_inc = RFIexcision.Compare_IB(puls_sb,puls_inc)
        
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
        inf.insert(0,'BEAM',str(beam))
        inf.insert(0,'SAP',str(sap))
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
    
    puls_inc_all = puls_inc_all.append(puls_inc,ignore_index=False)
      
  #Compares pulses in different beams
  puls = RFIexcision.Compare_Beams(puls)
  meta_data = meta_data.astype(str)
  
  puls = puls.append(puls_inc_all,ignore_index=False)
  puls_inc_all = 0
  
  best_puls = RFIexcision.best_pulses(puls,data)
  puls.sort(['SAP','BEAM','Pulse','Sigma'],ascending=[True,True,True,False],inplace=True)
  
  #Stores the table into a DataBase
  if not data.empty:
    store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'w')
    store.append(idL,data,data_columns=['Pulse'])
    store.append(idL+'_pulses',puls,data_columns=['Pulse'])
    store.append('meta_data',meta_data)
    if not best_puls.empty: store.append('best_pulses',best_puls)  #FORSE da togliere
    store.close()
    
  output(folder,idL,puls,best_puls,data,meta_data)
        
  return


def output(folder,idL,puls,best_puls,data,meta_data):

  if not data.empty:
    for sap in range(0,3):
      for beam in range(12,74):
        puls_plot = puls[(puls.SAP==sap)&(puls.BEAM==beam)]
        if not puls_plot.empty:
          name = 'SAP{}_BEAM{}'.format(sap,beam)
          os.makedirs('{}{}/sp/'.format(folder,idL)+name)
          #Pulse_min = puls_plot.Pulse.unique().min()
          astro = puls_plot[puls_plot.Pulse==0]
          rfi = puls_plot[puls_plot.Pulse==1]
          data_plot = data[data.Pulse.isin(astro.index)]
          meta_data_plot = meta_data[(meta_data.SAP==str(sap))&(meta_data.BEAM==str(beam))]
          best_puls_plot = best_puls[(best_puls.SAP==sap)&(best_puls.BEAM==beam)]
          if beam == 12: 
            LSPplot.plot(astro.iloc[30:],rfi,meta_data_plot,astro.iloc[:30],best_puls_plot,store='{}{}/sp/{}/{}_{}.png'.format(folder,idL,name,idL,name))
            LSPplot.sp(astro.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(0-9).png'.format(folder,idL,name))
            LSPplot.sp(astro.iloc[10:20],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(10-19).png'.format(folder,idL,name))
            LSPplot.sp(astro.iloc[20:30],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(20-29).png'.format(folder,idL,name))
            LSPplot.sp(best_puls_plot.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(0-9).png'.format(folder,idL,name))
            LSPplot.sp(best_puls_plot.iloc[10:20],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(10-19).png'.format(folder,idL,name))
            LSPplot.sp(best_puls_plot.iloc[20:30],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(20-29).png'.format(folder,idL,name))
          else:
            LSPplot.plot(astro.iloc[10:],rfi,meta_data_plot,astro.iloc[:10],best_puls_plot,store='{}{}/sp/{}/{}_{}.png'.format(folder,idL,name,idL,name))
            LSPplot.sp(astro.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/top_candidates.png'.format(folder,idL,name))
            LSPplot.sp(best_puls_plot,data,meta_data_plot,store='{}{}/sp/{}/best_pulses.png'.format(folder,idL,name))
    LSPplot.obs_top_candidates(puls[(puls.Pulse==0)&(puls.BEAM>12)].groupby(['SAP','BEAM'],sort=False).head(10),best_puls[best_puls.BEAM>12],store='{}{}/sp/top_candidates.png'.format(folder,idL)) 
    LSPplot.obs_top_candidates(puls[(puls.Pulse==0)&(puls.BEAM==12)].groupby('SAP',sort=False).head(30),best_puls[best_puls.BEAM==12],store='{}{}/sp/inc_top_candidates.png'.format(folder,idL),incoherent=True)
    
    #print puls[(puls.Pulse==0)].groupby(['SAP','BEAM'],sort=False).head(10)
    #print puls[(puls.Pulse==0)].groupby(['SAP','BEAM'],sort=False).max()
    
    best_puls['code'] = best_puls.index
    if not best_puls.empty:
      a = best_puls.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
      b = [val for sublist in a for val in sublist]
      best_puls.index = b
    best_puls.Duration *= 1000
    best_puls['void'] = ''
    best_puls.to_csv('{}{}/sp/best_pulses.inf'.format(folder,idL),sep='\t',float_format='%.2f',columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
    
    top_candidates = puls[puls.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
    top_candidates = top_candidates.append(puls[puls.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
    top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
    top_candidates['code'] = top_candidates.index
    if not top_candidates.empty:
      a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
      b = [val for sublist in a for val in sublist]
      top_candidates.index=b
    top_candidates.Duration *= 1000
    top_candidates['void'] = ''
    top_candidates.to_csv('{}{}/sp/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
        
  return
