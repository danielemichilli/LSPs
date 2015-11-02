import os
import numpy as np
import glob
import multiprocessing as mp


''' NEW

DM steps: 1

In standard pipeline:
- Applicare jumps removal
- Produrre 500 timeseries per ogni beam con -clip argument
- Chiamare beam_matrix
- Rimuovere timeseries
- Salvare matrice numpy
In launcher:





'''




def beam_matrix(DM_max,DM_step,timeseries_lenght):
  """
  Load the timesieries for a beam and transform the signal in SNR
  """
  
  #Load the timeseries in a matrix
  beam = np.zeros((DM_max/DM_step,timeseries_lenght),dtype=np.float32)
  
  #for file in os.listdir('.'):
    #if file.endswith(".dat"):
      #timeresies = np.fromfile(file,dtype=np.float32)
      #DM = os.path.splitext(os.path.basename(file))[0]
      #DM = DM.split('DM')[-1]
      #idx = np.int(np.round(float(DM)/DM_step))
      #beam[idx] = timeresies
      
  for idx,DM in enumerate(np.arange(0,DM_max,DM_step)):
    file = glob.glob('*_DM{}.dat'.format(DM))[0]
    beam[idx] = np.fromfile(file,dtype=np.float32)
  
  sharedBeam = mp.Array('f', beam, lock=False)
  
  #Convert signal in SNR
  N_cores = 24
  pool = mp.Pool()
  pool.map(time_fft, range(os.cpu_count()-1))
  pool.close()
  pool.join()
  
  time_fft(beam)
  
  np.save('../FRB_beam_matrix',beam)
      
  return
  
  
def time_fft(CPU):
  DM_range = int(DM_max)/os.cpu_count()


  for idx,DM in enumerate(np.arange(CPU*DM_range,(CPU+1)*DM_range,DM_step)):
    file = glob.glob('*_DM{}.dat'.format(DM))[0]
    beam[idx] = np.fromfile(file,dtype=np.float32)  
  
  chunk = sharedBeam[]
  downfacts = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
  chunklen = timeseries.shape[1]
  for idx in range(timeseries.shape[0]):
    ts = timeseries[idx]
    ts = signal.detrend(ts, type='constant')
    tmpchunk = ts.copy()
    tmpchunk.sort()
    stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
    stds *= 1.148
    ts /= stds
    
    down_series = np.zeros(ts.size)
    
    for downfact in downfacts:
      goodchunk = np.convolve(ts, np.ones(downfact) / np.sqrt(downfact), mode='same')
      down_series = np.max((goodchunk,down_series),axis=0)
      #store on a second array of ints the downfactor value for the highest SNR 
    
    ts = down_series
  
return np.unravel_index(np.argpartition(timeseries, -lim1, axis=None)[-lim1:],timeseries.shape)[0]
  
  
