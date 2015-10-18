import numpy as np



signal = #loaded from timeseries



def convolve(arr,max_length=1):
  arr = scipy.signal.detrend(arr, type='linear')
  chunklen = arr.size
  tmpchunk = arr.copy()
  tmpchunk.sort()
  stds = np.sqrt((tmpchunk[chunklen/40:-chunklen/40]**2.0).sum() / (0.95*chunklen))
  stds *= 1.148
  arr /= stds    

  if not isinstance(max_length,list):
    max_length = range(max_length)
  
  signal = np.zeros((chunklen,len(max_length)))
  
  for i,length in enumerate(max_length):
    signal[i] = np.convolve(arr, np.ones(length), mode='same') / np.sqrt(length)
  
  return np.sum(signal,axis=0)


