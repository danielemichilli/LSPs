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
import math

import Events
import Pulses
import RFIexcision
import LSPplot
from Parameters import *

def obs_events(folder,idL,load_events=False,conf=False):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------

  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')

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
  results = [lists_creation((folder,idL,sap,beam,store)) for (sap,beam) in folders]
  pool.close()
  pool.join()
  
  pulses = pd.concat(results)
  results = 0
  pulses.sort_index(inplace=True)

  store.append('pulses',pulses,data_columns=['Pulse'])
  store.close()
    
  #Compares pulses in different beams
  #RFIexcision.Compare_Beams(pulses)

  #Clean the pulses table
  #pulses = pulses[pulses.Pulse <= RFI_percent]
 
  #Clean the events table
  if pulses[pulses.Pulse==0].empty: logging.warning("Any reliable pulse detected!")
  else:  
    events = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'events',where=['Pulse==pulses.index.tolist()'])
    meta_data = pd.read_hdf('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'meta_data')
    
    #Produce the output
    pulses = pulses[pulses.Pulse <= 2]
    events = events[events.Pulse.isin(pulses.index)]
    output(folder,idL,pulses,events,meta_data)
  
  return



def lists_creation((folder,idL,sap,beam,store)):

  #Import the events
  events, meta_data = Events.Loader(folder,idL,sap,beam)
  store.append('meta_data',meta_data)
  
  pulses = pd.DataFrame()

  if not events.empty:
    try:
      #Correct for the time misalignment of events
      events.sort(['DM','Time'],inplace=True)
      events.Time = Events.TimeAlign(events.Time,events.DM)
      
      #Group the events
      events.sort(['DM','Time'],inplace=True)
      Events.Group(events)
      store.append('events',events,data_columns=['Pulse'])
      
      #Apply the thresholds to the events
      events = Events.Thresh(events)

      #Generate the pulses
      events = events[events.Pulse>=0]
      pulses = Pulses.Generator(events)
      
      #Apply RFI filters to the pulses
      pulses = pulses[pulses.Sigma >= 6.5]
      events = events[events.Pulse.isin(pulses.index)]
      RFIexcision.Pulse_Thresh(pulses,events)
      
      events = 0
      pulses = pulses[pulses.Pulse < RFI_percent]
      Pulses.Candidates(pulses)
      
    except:
      logging.warning("Some problem arised processing SAP "+str(sap)+" - BEAM "+str(beam)+", it will be discarded")
      return pd.DataFrame()

  return pulses



  


def output(folder,idL,pulses,events,meta_data):
  pulses.sort('Sigma',ascending=False,inplace=True)

  store = '{}{}/sp/beam_plots'.format(folder,idL)
  
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  
  for n in gb_puls.indices.iterkeys():
    name = 'SAP{}_BEAM{}'.format(n[0],n[1])
    os.makedirs('{}/{}'.format(store,name))
    
  pool = mp.Pool()
  pool.map(output_beams, [(pulses,events,meta_data,store,n) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()  
  
  #for n in gb_puls.indices.iterkeys():
    #LSPplot.plot((Group_Clean((gb_puls,gb_rfi,gb_md,gb_event),n),store,n))
  
  output_pointing(pulses,folder,idL)
  
  return



def output_beams((pulses,meta_data,store,(sap,beam))):
  
  top = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==0)]
  good = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==1)]
  
  if beam == 12:
    LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates(0-9).png'.format(folder,sap,beam),obs)
    LSPplot.sp_shape(top.iloc[10:20],events,'{}/SAP{}_BEAM{}/top_candidates(10-19).png'.format(folder,sap,beam),obs)
    LSPplot.sp_shape(top.iloc[20:30],events,'{}/SAP{}_BEAM{}/top_candidates(20-29).png'.format(folder,sap,beam),obs)
    LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
    
  else:
    LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates.png'.format(folder,sap,beam),obs)
    LSPplot.sp_plot(top,rfi,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
  
  return
  
  

def output_pointing(pulses,folder,idL):
  
  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index = b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}{}/sp/files/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',\
    columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],\
    header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  
  LSPplot.obs_top_candidates(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),\
                             store='{}{}/sp/files/inc_top_candidates.png'.format(folder,idL),incoherent=True)

  for sap in pulses.SAP.unique():
    LSPplot.obs_top_candidates(pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10),\
                               store='{}{}/sp/files/top_candidates(SAP{}).png'.format(folder,idL,sap))

  return








def alerts(pulses,folder,idL):
  store = '{}{}/sp/candidates'.format(folder,idL)
  os.makedirs('{}'.format(store))
  
  num = pulses.groupby(['SAP','BEAM'])['Pulse'].count().groupby(level=0).mean()
  span = 54500.*0.05/c(num,2)
  
  p=1/20 #check!
  
  
  def c(n,k):
    return math.factorial(n)/math.factorial(k)/math.factorial((n-k))

  #def diff(n,k):              
    #return 550*(100*c(n,k))**-1./(k-1)
  
  def p(n,k):              
    return c(n,k)/(9000./span)**(k-1)  #span: number of DMs
  
  def span(n,k):
    return 54500.*(10./222./c(n,k)**(1./(k-1)))  #span: number of DMs
  
  return
  
  
  
  
  
  

