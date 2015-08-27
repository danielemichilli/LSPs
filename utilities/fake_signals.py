spectra = read_filterbank(filename)
spectra[:] = 40



t0 = 1000 #bins
DM = 14.


import numpy as np
from Parameters import *

#filterbank file
freq = np.linspace(F_MIN,F_MAX,2592)
time = 4149 * DM * (np.power(freq,-2) - F_MAX**-2) / RES + t0

spectra[time.astype(np.int),np.arange(2591,-1,-1)] = 255





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

  spectra = np.memmap(filename, dtype=dtype, mode='r', offset=header_size, shape=(nspec, nchans))
  
  return spectra