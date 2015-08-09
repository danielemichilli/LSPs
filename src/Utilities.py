#import filterbank
#import os
#import numpy as np
#import sigproc
#import pyfits

from Parameters import *


def read_filterbank(filename,DM,bin_start):
  #Read data from filterbank file
  head = filterbank.read_header(filename)[0]

  nbits = head['nbits']
  nchans = head['nchans']

  dtype = filterbank.get_dtype(nbits)

  filfile = open(filename, 'rb')
  filfile.seek(0)
  paramname = ""
  while (paramname != 'HEADER_END'): paramname, val = sigproc.read_hdr_val(filfile)
  header_size = filfile.tell()
  filfile.close()

  file_size = os.stat(filename)[6]
  data_size = file_size - header_size
  bytes_per_spectrum= nchans * nbits / 8
  nspec = data_size / bytes_per_spectrum

  spectrum = np.memmap(filename, dtype=dtype, mode='r', offset=header_size, shape=(nspec, nchans))
  spectrum = np.fliplr(spectrum)

  bin_start -= DS_OFFSET
  bin_end = bin_start + DM2delay(DM) + DS_OFFSET
  
  spectrum = spectrum[bin_start:bin_end]
  ind = np.arange(0,2592,16)
  spectrum = np.delete(spectrum,ind,axis=1)
  
  return spectrum
 


def read_fits(filename,DM,bin_start,offset):
  bin_start = np.int(bin_start)
  fits = pyfits.open(filename,memmap=True)
  
  header = fits['SUBINT'].header
  
  N_channels = header['NCHAN']
  N_spectra = header['NSBLK']

  bin_start -= offset
  bin_end = bin_start + DM2delay(DM) + offset
  
  subint_start = bin_start/N_spectra
  subint_end = bin_end/N_spectra+1
  subint = fits['SUBINT'].data[subint_start:subint_end]['DATA']
  
  fits.close()

  subint = subint.reshape((subint_end-subint_start)*N_spectra,N_channels)
  
  subint = subint[bin_start%N_spectra:bin_start%N_spectra+bin_end]
  ind = np.arange(0,2592,16)
  subint = np.delete(subint,ind,axis=1)
  
  return subint, bin_start, bin_end
 
 

def read_header(filename):
  name, ext = os.path.splitext(filename)
  if ext == '.fits': 
    fits = pyfits.open(filename,memmap=True)
    header = fits['SUBINT'].header + fits['PRIMARY'].header
    fits.close()
    return header
  elif ext == '.fil':
    return filterbank.read_header(filename)[0]
  else: return None



def DM2delay(DM):
  delay = 4149. * (F_MIN**-2 - F_MAX**-2)
  return np.int(delay*DM/RES)


def delay4DM(DM):
  delay = 4149. * (F_MIN**-2 - F_MAX**-2)

  DM_steps = np.array((
          2.525,    5.055,    7.585,   10.115,   12.645,   15.175,
         17.705,   20.235,   22.765,   25.295,   27.825,   30.355,
         32.885,   35.415,   37.945,   40.475,   65.815,   91.115,
        116.415,  141.715,  242.965,  344.165,  445.365,  546.565))
  
  DM_n = (np.digitize([DM,], DM_steps) - 1)[0]
  
  if DM_n<9: DM_n = 0.253 * DM_n
  elif DM_n<15: DM_n = 0.253 * DM_n + 0.253
  elif DM_n<19: DM_n = 0.253 * 15 + 25.3 * delay * ( DM_n - 14 )
  else: DM_n = 0.253 * 15 + 25.3 * delay * 4 + 4 * 25.3 * delay * ( DM_n - 18 )
  
  return delay * DM / 2 + DM_n #s


def time2bin(time,DM):
  time -= delay4DM(DM)
  return time/RES


def span(n,k):
  return 54500.*(10./222./c(n,k)**(1./(k-1)))  #span: number of DMs
  
def c(n,k):
  return math.factorial(n)/math.factorial(k)/math.factorial((n-k))


def color_range(data):
  #Define the color range
  clean = data[data>0]
  min_element = clean.size/20
  max_element = clean.size*9/10
  vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
  vmax = np.partition(clean, max_element, axis=None)[max_element]
  return vmin,vmax