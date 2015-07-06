filename = 'L204720_SAP1_BEAM48.fil'
#pulsar_period = #s
#pulsar_DM = 


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import filterbank as fb
import os


#Open and elaborate the filterbank file
header, header_size = fb.read_header(filename)

f = open(filename, 'rb')
f.seek(header_size, os.SEEK_SET)
data = np.fromfile(f, dtype=np.uint8)
f.close()

#data = np.memmap(filename, dtype=np.uint8, offset=header_size)
#data = data.view(np.ma.MaskedArray)

nchans = header['nchans']
nspec = data.size/nchans

data = data.reshape((nspec, nchans))
data = np.fliplr(data)
ind = np.arange(0,2592,16)
data = np.delete(data,ind,axis=1)


#chunks = 16
#data = data.reshape(chunks,data.size/chunks)

#Define the 
tsamp = header['tsamp']
period = pulsar_period / tsamp
window = int(period)

offset = period-window

number_windows = int(data.size/window)

ind = []
for i in range(1,number_windows):
  offset *= i
  if offset >= 1: 
    ind.append(window*i)
    offset -= 1
    
np.delete(data,ind,axis=0)


