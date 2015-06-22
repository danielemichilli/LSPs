filename = 





import filterbank
import os
import numpy as np

head = filterbank.read_header(filename)[0]

nbits = head['nbits']
nchans = head['nchans']

dtype = filterbank.get_dtype(nbits)

filfile = open(filename, 'rb')
filfile.seek(0)
header_size = filfile.tell()
filfile.close()

file_size = os.stat(filename)[6]
data_size = file_size - header_size
bytes_per_spectrum= nchans * nbits / 8
nspec = data_size / bytes_per_spectrum

spectra = np.memmap(filename, dtype=dtype, mode='r', offset=header_size, shape=(nspec, nchans))

print spectra

