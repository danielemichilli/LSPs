import numpy as np
import glob
import multiprocessing as mp
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata


''' NEW

DM steps: 1

In standard pipeline:
- Applicare jumps removal
- Produrre 500 timeseries per ogni beam con -clip argument
- Chiamare beam_matrix
- Rimuovere timeseries
- Salvare matrici numpy
In LSPs:
- Chiamare FRBs_finder




'''

#datacube(61,500,100000)

Number_of_candidates = 10


def beam_matrix():
  '''
  Load the timesieries for a beam and transform the signal in SNR
  '''
  #Forse possibile detrend in DM ogni time bin?

  #Create and save the SNR datacube
  time_fft_parallel()
  
  #Create and save the signal datacube
  for idx,DM in enumerate(np.arange(0,500)):
    file = glob.glob('*_DM{}.dat'.format(DM))[0]
    if idx == 0: 
      beam0 = np.fromfile(file,dtype=np.float32)
      beam = np.zeros((500,beam0.size))
      beam[0] = beam0
    else:
      beam[idx] = np.fromfile(file,dtype=np.float32)
    
  np.save('../FRB_beam_matrix_signal',beam)  
  
  return
  
  
def time_fft_parallel():
  '''
  Parallelize the time_fft function
  '''
  CPUs = mp.cpu_count()
  pool = mp.Pool(CPUs)
  results = pool.map_async(time_fft, range(CPUs))
  pool.close()
  pool.join()
  
  beam = np.vstack(results)
  np.save('../FRB_beam_matrix_SNR',beam)
  
  return


def time_fft(CPU):
  '''
  Create the SNR datacube 
  '''
  CPUs = mp.cpu_count()
  DM_range = int(np.ceil(500/float()))  #Number of DMs in each CPU
  
  if DM_range * CPU > 500: timeseries = np.zeros((500%CPUs,36000))
  else: timeseries = np.zeros((DM_range,36000))
  downfacts = [2, 5, 10, 20, 50, 100, 200, 500, 1000]
  chunklen = timeseries.shape[1]
  
  #Produce the datacube
  for idx,DM in enumerate(np.arange(CPU*DM_range,np.clip((CPU+1)*DM_range,0,500))):
    file = glob.glob('*_DM{}.dat'.format(DM))[0]
    ts = np.fromfile(file,dtype=np.float32)[:36000]
    
    if idx == 0:
      timeseries = np.zeros((DM_range,ts.size))
      
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
    
    timeseries[idx] = down_series
  
  return timeseries




def dc_load():
  
  
  return datacube



def FRBs_finder():
  fill = 10  #Pixels per beam

  ra = np.array([399, 399, 498, 497, 399, 302, 301, 399, 499, 598, 595, 594, 496,
        399, 303, 206, 203, 201, 300, 399, 500, 600, 698, 695, 692, 690,
        592, 495, 399, 304, 207, 110, 107, 103, 101, 199, 299, 399, 502,
        602, 701, 800, 796, 792, 788, 785, 687, 590, 494, 399, 305, 209,
        112,  14,  11,   7,   3,   0,  98, 198, 298])

  dec = np.array([400, 500, 450, 350, 300, 350, 450, 600, 550, 500, 400, 300, 250,
        200, 250, 300, 400, 500, 550, 700, 650, 600, 550, 450, 350, 250,
        200, 150, 100, 150, 200, 250, 350, 450, 550, 600, 650, 800, 750,
        700, 650, 600, 500, 400, 300, 200, 150, 100,  50,   0,  50, 100,
        150, 200, 300, 400, 500, 600, 650, 700, 750])

  grid = np.linspace(0,800,8*fill)
  
  datacube = dc_load()
  
  #best_time_idxs = np.max(datacube,axis=(0,1)).argsort()[::-1]  #Assuming time is on axis 2
  
  cand_idx = np.zeros(Number_of_candidates) - 1 
  
  


def space_fft_parallel()
  pool = mp.Pool()
  results = pool.map(space_fft, range(mp.cpu_count()))
  pool.close()
  pool.join()
  
  beam = np.vstack(results)
  return beam
  
  
def space_fft(CPU):
  t_range = 36000 / (mp.cpu_count())
  
  beams = data[:,:,CPU*t_range:(CPU+1)*t_range]
  
  for time in times:
    #DM_idx = np.unravel_index(np.argmax(datacube[:,:,time_idx],datacube.shape)[1]
    
    ts = datacube[:,DM_idx,time_idx]
    ts = signal.detrend(ts, type='constant')
    tmpchunk = ts.copy()
    tmpchunk.sort()
    stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
    stds *= 1.148
    ts /= stds

    beams_map = griddata(ra, dec, ts, grid, grid, interp='nn')
    beams_map[beams_map.mask==True] = 0
    conv = signal.convolve2d(beams_map,np.ones((fill,fill))/fill,mode='same')  #Check the statistics!
    
    
    #Condizioni per salvare il bin


    if cand_idx[cand_idx>=0].size == 0: break
  
  
#Create the convolved matrix
beams_map = np.zeros((9,9))
kern1 =  = np.ones((2,3))
kern2 =  = np.ones((3,2))
conv1 = signal.convolve2d(beams_map,kern1,mode='full')
conv2 = signal.convolve2d(beams_map,kern2,mode='full')

conv1[1:-1,1] /= 2
conv1[1:-1,-2] /= 2
conv1[0,2:-2:2] /= 2
conv1[-1,2:-2:2] /= 2
conv1[1:-1,2:-2] /= 3

conv2[1,1:-1] /= 2
conv2[-2,1:-1] /= 2
conv2[2:-2:2,0] /= 2
conv2[2:-2:2,-1] /= 2
conv2[2:-2,1:-1] /= 3

conv = #Bisogna prendere il massimo per ogni beam