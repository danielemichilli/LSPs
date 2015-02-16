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

import Events
import Pulses
import RFIexcision
import LSPplot
from Parameters import *


def Initialize():  #TOGLIERE!
  #Creates the tables in memory
  meta_data = pd.DataFrame(columns=['SAP','BEAM','File','Telescope','Instrument','RA','DEC','Epoch'])
  events = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse'])
  pulses = pd.DataFrame(columns=['SAP','BEAM','DM','Sigma','Time','Duration','Pulse','dDM','dTime','DM_c','Time_c','N_events'])
  
  pulses = pulses.astype(np.float32)
  pulses.SAP = pulses.SAP.astype(np.uint8)
  pulses.BEAM = pulses.BEAM.astype(np.uint8)
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses.N_events = pulses.N_events.astype(np.int16)
  
  meta_data = meta_data.astype(str)
  meta_data.SAP = meta_data.SAP.astype(np.uint8)
  meta_data.BEAM = meta_data.BEAM.astype(np.uint8)
  
  return events,pulses,meta_data



def obs_events(folder,idL):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------
  
    
  
  #Initialize the tables
  #events,pulses,meta_data = Initialize()
  
  #Import the events
  pool = mp.Pool(mp.cpu_count()-1)
  results = pool.map(Events.Loader, [(folder,idL,beam) for beam in range(12,74)])
  pool.close()
  pool.join()
  
  results_events = [x[0] for x in results]
  results_meta_data = [x[1] for x in results]
  results = 0
  
  events = pd.concat(results_events)
  meta_data = pd.concat(results_meta_data)
  results_events = 0
  results_meta_data = 0
  
  #events, meta_data = Events.Loader((folder,idL,12))
  
  events['Pulse'] = 0
  
  events = events.astype(np.float32)
  events.SAP = events.SAP.astype(np.uint8)
  events.BEAM = events.BEAM.astype(np.uint8)
  events.Pulse = events.Pulse.astype(np.int32)
  
  #Apply the thresholds to the events
  events = Events.Thresh(events)

    
  #Group the events
  Events.Group(events)
  
  events = events[events.Pulse>0]
    
  #Generate the pulses table  
  rows_core = events.shape[0] / (mp.cpu_count()+2)
  
  gb = events.groupby('BEAM',sort=False)
  
  pool = mp.Pool(mp.cpu_count()-1)
  results = pool.map(Pulses.Generator, [gb.get_group(n) for n in gb.indices.iterkeys()])
  pool.close()
  pool.join()
  #output = [p.get() for p in results] 
  
  pulses = pd.concat(results)
  results = 0
    
  #pulses = Pulses.Generator(pulses,events)
  
  
  #Correct for the time misalignment of events
  Pulses.TimeAlign(events.Time,events.DM)
  
  #Store the events
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'w')
  store.append(idL,events,data_columns=['Pulse'])
  store.close() 
  
  #Clean the events table
  events = events.loc[events.Pulse.isin(pulses.index)]
  
  #Compares pulses in different beams
  RFIexcision.Compare_Beams(pulses)

  #Clean the pulses table
  pulses = pulses[pulses.Pulse <= RFI_percent]







  
  def ALERT(sap,beam,dm,limit,counts,s_max,s_mean):
    file = open('{}{}/sp/ALERTS'.format(folder,idL),'a+')
    if not file.readlines():
      file.write('Pulsar candidates\n\n')
      file.write('SAP\tBEAM\tDM\tLimit\tCounts\tS_Max\tS_Mean\n')
    file.write('{0}\t{1}\t{2:.2f}\t{3:.2f}\t{4}\t{5:.2f}\t{6:.2f}\n'.format(sap,beam,dm,limit,counts,s_max,s_mean))
    file.close()
    return
  
  
  
  
  data = pulses.loc[pulses.BEAM>12,['SAP','BEAM','DM','Sigma','Pulse']]
  data.DM = data.DM.round(decimals=1)
  noise_level = data.groupby(['SAP','BEAM','DM'],sort=False).Pulse.count()
  noise_level = noise_level.mean(level=['SAP','DM']) + 5.* noise_level.std(level=['SAP','DM'])   #5 sigmas tollerance
  noise_level.dropna(inplace=True)
    
  beams = data.groupby(['SAP','BEAM'],sort=False)
  for ind,beam in beams:
    counts = beam.groupby(['SAP','DM'],sort=False).Pulse.count()
    counts = counts[(counts>noise_level.loc[counts.index])&(counts>5)]
    if not counts.empty:
      for i in counts.index.get_level_values('DM'):
        Sigma = beams.get_group((ind[0],ind[1]))
        Sigma = Sigma.Sigma[Sigma.DM==i]
        ALERT(ind[0],ind[1],i,noise_level.loc[ind[0],i],counts.loc[ind[0],i],Sigma.max(),Sigma.mean())
      



        
        
  
  # ALERT su incoherent beams rispetto a cosa? media dei tre? media di piu' osservazioni?
  
  #noise_level = pulses.loc[pulses.BEAM>12].groupby(['SAP','BEAM','DM'],sort=False).count().Pulse
  #noise_level = noise_level.mean(level='DM')# + 3.* noise_level.std(level='DM')   #3 sigmas tollerance
  #noise_level.dropna(inplace=True)
  
  #beams = pulses.loc[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False)
  #for ind,beam in beams:
    #counts = beam.groupby('DM',sort=False).Pulse.count()
    #counts = counts[counts>noise_level.loc[counts.index]]
    #if not counts.empty:
      #ALERT(beam.SAP.iloc[0],beam.BEAM.iloc[0],counts.index.values)




 
  
  #Generate best_puls table
  best_puls = RFIexcision.best_pulses(pulses,events)

  #Store the pulses
  store = pd.HDFStore('{}{}/sp/SinglePulses.hdf5'.format(folder,idL),'a')
  store.append(idL+'_pulses',pulses,data_columns=['Pulse'])
  store.append('meta_data',meta_data)
  if not best_puls.empty: store.append('best_pulses',best_puls)
  store.close()
  
  #Produce the output
  output(folder,idL,pulses,best_puls,events,meta_data)
  
  return


def output(folder,idL,pulses,best_puls,events,meta_data):
  
  store = '{}{}'.format(folder,idL)
  
  pulses.sort('Sigma',ascending=False,inplace=True)
  best_puls.sort('Sigma',ascending=False,inplace=True)
  
  gb_rfi = pulses.loc[pulses.Pulse==1].groupby(['SAP','BEAM'],sort=False)
  
  pulses = pulses.loc[pulses.Pulse==0]
    
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  
  gb_best = best_puls.groupby(['SAP','BEAM'],sort=False)

  gb_event = events.loc[(events.Pulse.isin(gb_puls.head(30).index))|(events.Pulse.isin(gb_best.head(30).index))].groupby(['SAP','BEAM'],sort=False)
  events = 0  
  
  gb_md = meta_data.groupby(['SAP','BEAM'],sort=False)
  meta_data = 0
  
  for n in gb_puls.indices.iterkeys():
    name = 'SAP{}_BEAM{}'.format(n[0],n[1])
    os.makedirs('{}{}/sp/'.format(folder,idL)+name)
  
  
  
  def Group_Clean(gb_best,n):
    try:
      return gb_best.get_group(n)
    except KeyError:
      return pd.DataFrame()
  
  
  #METTERE che n salta best_pulses se non sono presenti  
  pool = mp.Pool(mp.cpu_count()-1)
  pool.map(LSPplot.plot, [(gb_puls.get_group(n),gb_rfi.get_group(n),gb_md.get_group(n),Group_Clean(gb_best,n),gb_event.get_group(n),store) for n in gb_puls.indices.iterkeys()])
  pool.close()
  pool.join()
  
  
  best_puls['code'] = best_puls.index
  if not best_puls.empty:
    a = best_puls.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    best_puls.index = b
  best_puls.Duration *= 1000
  best_puls['void'] = ''
  best_puls.to_csv('{}{}/sp/best_pulses.inf'.format(folder,idL),sep='\t',float_format='%.2f',columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index=b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}{}/sp/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  
  LSPplot.obs_top_candidates(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),best_puls[best_puls.BEAM==12],\
                             store='{}{}/sp/inc_top_candidates.png'.format(folder,idL),incoherent=True)

  for sap in pulses.SAP.unique():
    LSPplot.obs_top_candidates(pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10),best_puls[(best_puls.SAP==sap)&(best_puls.BEAM>12)],\
                               store='{}{}/sp/top_candidates(SAP{}).png'.format(folder,idL,sap)) 

  pulses = 0

  
  #LSPplot.plot(astro.iloc[10:],rfi,meta_data_plot,astro.iloc[:10],best_puls_plot,store='{}{}/sp/{}/{}_{}.png'.format(folder,idL,name,idL,name))
  
  #for n in gb.indices.iterkeys():
    #puls_plot = gb.get_group(n)
    #if not puls_plot.empty:
      #name = 'SAP{}_BEAM{}'.format(n[0],n[1])
      #os.makedirs('{}{}/sp/'.format(folder,idL)+name)
      ##Pulse_min = puls_plot.Pulse.unique().min()
      #astro = puls_plot[puls_plot.Pulse==0]
      #rfi = puls_plot[puls_plot.Pulse==1]
      #data_plot = data[data.Pulse.isin(astro.index)]
      #meta_data_plot = meta_data[(meta_data.SAP==str(n[0]))&(meta_data.BEAM==str(n[1]))]
      #best_puls_plot = best_puls[(best_puls.SAP==n[1])&(best_puls.BEAM==n[1])]
      #if n[1] == 12: 
        #LSPplot.plot(astro.iloc[30:],rfi,meta_data_plot,astro.iloc[:30],best_puls_plot,store='{}{}/sp/{}/{}_{}.png'.format(folder,idL,name,idL,name))
        #LSPplot.sp(astro.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(0-9).png'.format(folder,idL,name))
        #LSPplot.sp(astro.iloc[10:20],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(10-19).png'.format(folder,idL,name))
        #LSPplot.sp(astro.iloc[20:30],data,meta_data_plot,store='{}{}/sp/{}/top_candidates(20-29).png'.format(folder,idL,name))
        #LSPplot.sp(best_puls_plot.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(0-9).png'.format(folder,idL,name))
        #LSPplot.sp(best_puls_plot.iloc[10:20],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(10-19).png'.format(folder,idL,name))
        #LSPplot.sp(best_puls_plot.iloc[20:30],data,meta_data_plot,store='{}{}/sp/{}/best_pulses(20-29).png'.format(folder,idL,name))
      #else:
        #LSPplot.plot(astro.iloc[10:],rfi,meta_data_plot,astro.iloc[:10],best_puls_plot,store='{}{}/sp/{}/{}_{}.png'.format(folder,idL,name,idL,name))
        #LSPplot.sp(astro.iloc[:10],data,meta_data_plot,store='{}{}/sp/{}/top_candidates.png'.format(folder,idL,name))
        #LSPplot.sp(best_puls_plot,data,meta_data_plot,store='{}{}/sp/{}/best_pulses.png'.format(folder,idL,name))

   
    #print puls[(puls.Pulse==0)].groupby(['SAP','BEAM'],sort=False).head(10)
    #print puls[(puls.Pulse==0)].groupby(['SAP','BEAM'],sort=False).max()
    
  return
