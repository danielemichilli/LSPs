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
import fnmatch
import matplotlib as mpl
mpl.use('Agg')
import multiprocessing as mp
import matplotlib.pyplot as plt
import logging
import re

import Events
import Pulses
import RFIexcision
import LSPplot
import Candidates
import Internet
from Parameters import *
from Paths import *


def obs_events(args, debug=False):
  #----------------------------------------------------------
  # Creates the clean table for one observation and stores it
  #----------------------------------------------------------
  if args.conf: inc = 0
  else: inc = 12
  
  def file_list(folder, idL):
    file_list = []
    file_names = []
    for root, dirnames, filenames in os.walk(os.path.join(folder+idL)):
      for filename in fnmatch.filter(filenames, '*singlepulse.tgz'):
        if filename in file_names:
          logging.warning("ATTENTION!: Two sp archives with the same filename foud: {}. Only the latter will be processed".format(filename))
          idx = file_names.index(filename)
          file_list[idx] = root+'/'+filename
          continue
        else:
          file_list.append(root+'/'+filename)
          file_names.append(filename)
    
    return file_list
  
  file_list = file_list(args.folder,args.idL)
    
  if debug:
    CPUs = mp.cpu_count()
    dirs_range = int(np.ceil(len(file_list)/float(CPUs)))
    results = [lists_creation(args.idL,file_list[i*dirs_range:(i+1)*dirs_range],args.folder) for i in range(CPUs) if len(file_list[i*dirs_range:(i+1)*dirs_range]) > 0]
    pulses = pd.concat(results)
  else: pulses = pulses_parallel(args.idL,file_list,args.folder)
  
  def merge_temp_databases(idL,store,file):
    store.append('events',pd.read_hdf('{}/{}'.format(TMP_FOLDER.format(idL),file),'events'),data_columns=['Pulse','SAP','BEAM','DM','Time'])
    meta_data = pd.read_hdf('{}/{}'.format(TMP_FOLDER.format(idL),file),'meta_data')
    meta_data.reset_index(inplace=True,drop=True)
    store.append('meta_data',meta_data)    
    os.remove('{}/{}'.format(TMP_FOLDER.format(idL),file))
    
  store = pd.HDFStore('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(args.idL)),'w')
  features_list = ''
  for file in os.listdir(TMP_FOLDER.format(args.idL)):
    if file.endswith('.arff_tmp'):
      with open(os.path.join(TMP_FOLDER.format(args.idL),file), 'r') as f:
        line = f.readline()
        idx = len(line.split(',')) - 1
        break
  
  for i in range(idx): features_list += '@attribute Feature{} numeric\n'.format(i)
  header = """@relation Training_v3
  {}
  @attribute class {{0,1}}
  @data
  """.format(features_list[:-1])
  thresholds = open('{}/thresholds.arff'.format(TMP_FOLDER.format(args.idL)), 'w')
  thresholds.write(header)
  for file in os.listdir(TMP_FOLDER.format(args.idL)):
    if file.endswith('.tmp'):
      merge_temp_databases(args.idL,store,file)
    if file.endswith('.arff_tmp'):
      with open(os.path.join(TMP_FOLDER.format(args.idL),file), 'r') as f:
        thresholds.write(f.read())
      os.remove(os.path.join(TMP_FOLDER.format(args.idL),file))
  thresholds.close()
  store.close()
  
  #Select positive pulses
  ML_predict = os.path.join(TMP_FOLDER.format(args.idL), 'ML_predict.txt')  
  pulses = select_real_pulses(pulses,'{}/thresholds'.format(TMP_FOLDER.format(args.idL)), ML_predict)
  
  
  if pulses.empty: 
    logging.warning("No pulse detected!")
    return
  pulses.sort_index(inplace=True)



  #TO BE TESTED

  '''
  # Nuove idee: 
  # filtro 1 - rimuove pulses in cui, escludedno i beam vicini, ci sono piu di 1/10 di pulses in altri beam in dt, dDM
  # filtro 2 - rimuove pulses in cui ci sono piu' di 1/10 di beams aventi la meta' del SNR del beam con il massimo segnale entro dt, dDM
  # filtro 3 - raggruppa i pulses vicini in t,DM e applica statistiche simili agli eventi con i pulses. nel raggruppamento, se ci sono due pulses entro dt,dDM si raggruppa il piu' vicino
  
  
  #Remove time-spans affectd by RFI
  pulses.sort_index(inplace=True)
  pulses.Pulse.loc[RFIexcision.time_span(pulses[pulses.BEAM == inc])] += 1
  pulses.Pulse.loc[RFIexcision.time_span(pulses[pulses.BEAM > inc])] += 1
  pulses = pulses[pulses.Pulse <= RFI_percent]

  #Compare different beams
  pulses.Pulse.loc[RFIexcision.Compare_Beams(pulses[pulses.BEAM > inc].copy())] += 1
  pulses = pulses[pulses.Pulse <= RFI_percent]

  #Remove pulses appearing in too many beams
  pulses.Pulse += pulses.apply(lambda x: RFIexcision.puls_beams_select(x,pulses),axis=1).astype(np.int8)
  pulses = pulses[pulses.Pulse <= RFI_percent]
  '''
  
  pulses.Candidate = pulses.Candidate.astype(np.int32)
  cands = Candidates.candidates(pulses,args.idL)
  
  store = pd.HDFStore('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(args.idL)),'a')
  store.append('pulses',pulses)
  if not cands.empty:
    cands.sort(['Sigma','Rank'],ascending=[0,1],inplace=True)
    store.append('candidates',cands)
    #cands = cands[cands.main_cand==0].head(30)
  store.close()
      
  if cands.empty: logging.warning("Any reliable candidate detected!")
  else:
    #Produce the output
    meta_data = pd.read_hdf('{}/sp/SinglePulses.hdf5'.format(WRK_FOLDER.format(args.idL)),'meta_data')
    LSPplot.output(args.idL, pulses, meta_data, cands, inc=inc)    
    
    #Store the best candidates online
    try: Internet.upload(cands,args.idL,'{}/sp/candidates/.'.format(WRK_FOLDER.format(args.idL)),meta_data)
    except: 
      logging.exception("ATTENTION!\n\nConnession problem, update candidates in a second moment\n\n")
      with open('{}/{}/SP_ERROR.txt'.format(args.folder,args.args.idL),'a') as f:
        f.write("Connession problem \nConsider to run Upload.py script\n")

  return


def pulses_parallel(idL,dirs,folder):    
  #Create events, meta_data and pulses lists
  CPUs = mp.cpu_count()
  dirs_range = int(np.ceil(len(dirs)/float(CPUs)))
  
  pool = mp.Pool(CPUs)
  results = [pool.apply_async(lists_creation, args=(idL,dirs[i*dirs_range:(i+1)*dirs_range],folder)) for i in range(CPUs) if len(dirs[i*dirs_range:(i+1)*dirs_range]) > 0]
  pool.close()
  pool.join()
  return pd.concat([p.get() for p in results])


def lists_creation(idL,dirs,folder):
  result = pd.DataFrame()
  
  for directory in dirs:
    def coord_from_path(directory):
      elements = os.path.basename(directory).split('_')
      sap = int(re.findall('\d', elements[1])[0])
      beam = int(re.findall('\d+', elements[2])[0])
      return sap, beam
    sap, beam = coord_from_path(directory)

    #try:
      #pulses = load_pulses(idL, directory, sap, beam)
      
    #except:
      #logging.warning("Some problem arised processing SAP "+str(sap)+" - BEAM "+str(beam)+", it will be discarded")
      #with open('{}/{}/SP_ERROR.txt'.format(folder,idL),'a') as f:
        #f.write("SAP {} - BEAM {} not processed due to some unknown error\n".format(sap, beam))
      #pulses = pd.DataFrame()
      
    result = pd.concat((load_pulses(idL, directory, sap, beam),result))
  
  return result


def load_pulses(idL, directory, sap, beam):
  #Load the events
  events, meta_data = Events.Loader(directory,sap,beam)
  pulses = pd.DataFrame()
  
  if events.empty: return pd.DataFrame()

  #Store the events        
  events.to_hdf('{}/SAP{}_BEAM{}.tmp'.format(TMP_FOLDER.format(idL),sap,beam),'events',mode='w')
  meta_data.to_hdf('{}/SAP{}_BEAM{}.tmp'.format(TMP_FOLDER.format(idL),sap,beam),'meta_data',mode='a')

  return pulses_from_events(events, idL, sap, beam)
  

def pulses_from_events(events, idL, sap, beam):
  #Correct for the time misalignment of events
  events.sort(['DM','Time'],inplace=True)  #Needed by TimeAlign
  events.Time = Events.TimeAlign(events.Time.copy(),events.DM)

  #Apply the thresholds to the events
  events = Events.Thresh(events)
  
  #Group the events
  events.sort(['DM','Time'],inplace=True) #Needed by Group
  Events.Group(events)
  events = events[events.Pulse>=0]

  #Generate the pulses
  pulses = Pulses.Generator(events)
  pulses = pulses[pulses.Sigma >= 6.5]

  #Set a maximum amout of pulses to prevent bad observations to block the pipeline
  pulses.sort('Sigma',ascending=False,inplace=True)
  pulses = pulses.iloc[:3e4]
  events = events[events.Pulse.isin(pulses.index)]
  
  #Apply RFI filters to the pulses
  if not pulses.empty:
    arff_basename = '{}/thresholds_{}_{}.arff_tmp'.format(TMP_FOLDER.format(idL), sap, beam)
    RFIexcision.filters(pulses, events, arff_basename, header=False)
 
  ##Remove weaker pulses within a temporal window
  #def simultaneous(p):                            
    #puls = pulses.Pulse[np.abs(pulses.Time-p.Time) < 0.02]
    #if puls.shape[0] == 1: return 0
    #if p.name == puls.index[0]: return 0
    #else: return 1
  #pulses.Pulse += pulses.apply(lambda x: simultaneous(x), axis=1)
  
  
  #A set of known pulses is necessary
  #Apply multimoment analysis to the pulses
  #events = events[events.Pulse.isin(pulses.index)]
  #RFIexcision.multimoment(pulses,idL)
  #pulses = pulses[pulses.Pulse <= RFI_percent]
  #events = events[events.Pulse.isin(pulses.index)]

  return pulses

