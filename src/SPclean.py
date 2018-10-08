#############################
# Single Pulse cleaner
# Written by Daniele Michilli
#############################

import os
import fnmatch
import multiprocessing as mp
import re

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import Events
import Pulses
import RFIexcision
import LSPplot
import Candidates
import Internet
from Paths import *


import time



def main(args):

  time0 = time.time()

  if args.conf: inc = 0
  else: inc = 12
  
  def file_list(folder, id_obs):
    file_list = []
    file_names = []
    for root, dirnames, filenames in os.walk(os.path.join(folder, id_obs)):
      for filename in fnmatch.filter(filenames, '*singlepulse.tgz'):
        if filename in file_names:
          log_err("Two sp archives with the same filename foud: {}. Only the latter will be processed".format(filename), id_obs)
          idx = file_names.index(filename)
          file_list[idx] = os.path.join(root,filename)
          continue
        else:
          file_list.append(os.path.join(root,filename))
          file_names.append(filename)
    return file_list
  
  file_list = file_list(args.folder,args.id_obs)
    

  print "time 0:", time.time() - time0
  time0 = time.time()



  if args.debug:
    CPUs = mp.cpu_count()
    dirs_range = int(np.ceil(len(file_list)/float(CPUs)))
    results = [lists_creation(args.id_obs,file_list[i*dirs_range:(i+1)*dirs_range],args.folder) for i in range(CPUs) if len(file_list[i*dirs_range:(i+1)*dirs_range]) > 0]
    pulses = pd.concat(results)
  else: pulses = pulses_parallel(args.id_obs,file_list,args.folder)


  print "time 1:", time.time() - time0
  time0 = time.time()

  
  def merge_temp_databases(id_obs,store,file):
    store.append('events',pd.read_hdf(os.path.join(TMP_FOLDER.format(id_obs),file),'events'),data_columns=['Pulse','SAP','BEAM','DM','Time'])
    meta_data = pd.read_hdf(os.path.join(TMP_FOLDER.format(id_obs),file),'meta_data')
    meta_data.reset_index(inplace=True,drop=True)
    meta_data['version'] = args.vers
    store.append('meta_data',meta_data)
    os.remove(os.path.join(TMP_FOLDER.format(id_obs),file))
    
  store = pd.HDFStore(os.path.join(WRK_FOLDER.format(args.id_obs),'sp/SinglePulses.hdf5'),'w')
  features_list = ''
  for file in os.listdir(TMP_FOLDER.format(args.id_obs)):
    if file.endswith('.arff_tmp'):
      with open(os.path.join(TMP_FOLDER.format(args.id_obs),file), 'r') as f:
        line = f.readline()
        idx = len(line.split(',')) - 1
        break


  print "time 2:", time.time() - time0
  time0 = time.time()



  for i in range(idx): features_list += '@attribute Feature{} numeric\n'.format(i)
  header = """@relation Training_v3
{}
@attribute class {{0,1}}
@data
  """.format(features_list[:-1])
  thresholds = open(os.path.join(TMP_FOLDER.format(args.id_obs),'thresholds.arff'), 'w')
  thresholds.write(header)
  for file in os.listdir(TMP_FOLDER.format(args.id_obs)):
    if file.endswith('.tmp'):
      merge_temp_databases(args.id_obs,store,file)
    if file.endswith('.arff_tmp'):
      with open(os.path.join(TMP_FOLDER.format(args.id_obs),file), 'r') as f:
        thresholds.write(f.read())
      os.remove(os.path.join(TMP_FOLDER.format(args.id_obs),file))
  thresholds.close()
  store.close()
  
  #Select positive pulses
  print "Total pulses produced: {}".format(pulses.shape[0])
  ML_predict = os.path.join(TMP_FOLDER.format(args.id_obs), 'ML_predict.txt')  
  pulses = RFIexcision.select_real_pulses(pulses,os.path.join(TMP_FOLDER.format(args.id_obs),'thresholds'), ML_predict)
  print "Pulses positively classified: {}".format(pulses.shape[0])


  print "time 3:", time.time() - time0
  time0 = time.time()



  pulses = RFIexcision.beam_comparison(pulses,database=os.path.join(WRK_FOLDER.format(args.id_obs),'sp/SinglePulses.hdf5'), inc=inc)
  print "Pulses after beam comparison: {}".format(pulses.shape[0])


  print "time 4:", time.time() - time0
  time0 = time.time()


  
  if pulses.empty: 
    print "No pulse detected!"
    log("No pulse detected!", args.id_obs)
    return
  pulses.sort_index(inplace=True)

  pulses.Candidate = pulses.Candidate.astype(np.int32)
  cands = Candidates.candidates(pulses,args.id_obs)
  
  store = pd.HDFStore(os.path.join(WRK_FOLDER.format(args.id_obs),'sp/SinglePulses.hdf5'),'a')
  store.append('pulses',pulses)
  if not cands.empty:
    store.append('candidates',cands)
  store.close()
      
  if cands.empty:
    print "Any reliable candidate detected!"
    log("Any reliable candidate detected!", args.id_obs)
    return

  cands = cands[cands.main_cand == 0]
  cands.sort_values('Sigma', inplace=True, ascending=False)
  cands = cands.groupby('BEAM').head(10)
  cands = cands.head(50)
  #best_cands = cands[cands.N_pulses==1].groupby('BEAM').head(2).groupby('SAP').head(4)  #Select brightest unique candidates, 2 per BEAM and 4 per SAP
  #best_cands = best_cands.append(cands[cands.N_pulses>1].groupby('BEAM').head(2).groupby('SAP').head(6))  #Select brightest unique candidates, 2 per BEAM and 6 per SAP
  cands = cands[ ((cands.N_pulses == 1) & (cands.Sigma>10.)) | ((cands.N_pulses > 1) & (cands.Sigma>16.)) ]
  cands.sort_values('Sigma', inplace=True, ascending=False)



  print "time 5:", time.time() - time0
  time0 = time.time()



  #pulses = pulses[pulses.Candidate.isin(cands.index)]
  #Produce the output
  meta_data = pd.read_hdf(os.path.join(WRK_FOLDER.format(args.id_obs),'sp/SinglePulses.hdf5'),'meta_data')
  LSPplot.output(args.id_obs, pulses, meta_data, cands, os.path.join(WRK_FOLDER.format(args.id_obs),'sp/SinglePulses.hdf5'), inc=inc)


  print "time 6:", time.time() - time0
  time0 = time.time()


  
  #Store the best candidates online
  try: Internet.upload(cands,args.id_obs,os.path.join(WRK_FOLDER.format(args.id_obs),'sp/candidates/.'),meta_data,pulses)
  except: log_err("Connession problem \nConsider to run Upload.py script\n", args.id_obs)



  print "time 7:", time.time() - time0
  time0 = time.time()



  return


def pulses_parallel(id_obs,dirs,folder):    
  #Create events, meta_data and pulses lists
  CPUs = mp.cpu_count()

  print "CPUs:", CPUs

  dirs_range = int(np.ceil(len(dirs)/float(CPUs)))
  
  pool = mp.Pool(CPUs)
  results = [pool.apply_async(lists_creation, args=(id_obs,dirs[i*dirs_range:(i+1)*dirs_range],folder)) for i in range(CPUs) if len(dirs[i*dirs_range:(i+1)*dirs_range]) > 0]
  pool.close()
  pool.join()
  return pd.concat([p.get() for p in results])


def lists_creation(id_obs,dirs,folder):
  result = pd.DataFrame()
  
  for directory in dirs:
    def coord_from_path(directory):
      elements = os.path.basename(directory).split('_')
      sap = int(re.findall('\d', elements[1])[0])
      beam = int(re.findall('\d+', elements[2])[0])
      return sap, beam
    sap, beam = coord_from_path(directory)
      
    result = pd.concat((pulses_from_events(id_obs, directory, sap, beam),result))
  
  return result


def pulses_from_events(id_obs, directory, sap, beam):
  #Load the events
  events, meta_data = Events.Loader(directory,sap,beam)
  pulses = pd.DataFrame()
  
  if events.empty: return pd.DataFrame()
  
  #Correct for the time misalignment of events
  events['Time_org'] = events.Time.copy()
  events.sort_values(['DM','Time'],inplace=True)  #Needed by TimeAlign
  events.Time = Events.TimeAlign(events.Time.copy(),events.DM)

  #Apply the thresholds to the events
  events = Events.Thresh(events)
  
  #Group the events
  events.sort_values(['DM','Time'],inplace=True) #Needed by Group
  Events.Group(events)
  events = events[events.Pulse>=0]

  #Store the events
  events.sort_values(['DM','Time'],inplace=True)
  events.to_hdf(os.path.join(TMP_FOLDER.format(id_obs),'SAP{}_BEAM{}.tmp'.format(sap,beam)),'events',mode='w')  #Deve andare prima!
  meta_data.to_hdf(os.path.join(TMP_FOLDER.format(id_obs),'SAP{}_BEAM{}.tmp'.format(sap,beam)),'meta_data',mode='a')

  #Generate the pulses
  pulses = Pulses.Generator(events)
  pulses = pulses[pulses.Sigma >= 6.5]
  pulses = pulses[pulses.DM >= 3.]

  #Set a maximum amout of pulses to prevent bad observations to block the pipeline
  pulses.sort_values('Sigma',ascending=False,inplace=True)
  pulses = pulses.iloc[:int(3e4)]
  events = events[events.Pulse.isin(pulses.index)]
  
  #Apply RFI filters to the pulses
  if not pulses.empty:
    arff_basename = os.path.join(TMP_FOLDER.format(id_obs),'thresholds_{}_{}.arff_tmp'.format(sap, beam))
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
  #RFIexcision.multimoment(pulses,id_obs)
  #pulses = pulses[pulses.Pulse <= RFI_percent]
  #events = events[events.Pulse.isin(pulses.index)]

  return pulses


def log(text, id_obs):
  file_out='{}/sp/log.txt'.format(WRK_FOLDER.format(id_obs))
  with open(file_out, 'a') as f:
    f.write(text)
    f.write('\n')
  return


def log_err(text, id_obs):
  file_out='{}/{}/SP_ERROR.txt'.format(OBS_FOLDER, id_obs)
  with open(file_out, 'a') as f:
    f.write(text)
    f.write('\n')
  return

