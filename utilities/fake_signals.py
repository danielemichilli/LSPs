import numpy as np
from Parameters import *

spectra = read_filterbank(filename)

t0 = 1000 #bins
DM = 12.5
duration = 20 #bins

#filterbank file
freq = np.linspace(F_MIN,F_MAX,2592)
time = (4149 * DM * (np.power(freq,-2) - F_MAX**-2) / RES).round().astype(np.int) + t0

for i in range(duration):
  spectra[time+i,np.arange(2591,-1,-1)] = np.random.randint(0,255,2592)

ind = np.arange(15,2592,16)
spectra[:10000,ind] = 0




def read_filterbank(filename):
  import filterbank
  import os
  import numpy as np
  import sigproc
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
  spectra = np.memmap(filename, dtype=dtype, mode='r+', offset=header_size, shape=(nspec, nchans))
  return spectra

