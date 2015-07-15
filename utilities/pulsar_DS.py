filename = 'L204720_SAP1_BEAM48.fil'
pulsar_period = 0.714519699726 #s
pulsar_DM = 26.7641


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

#chunks = 16
#data = data.reshape(chunks,data.size/chunks)


#Define the process parameters
tsamp = header['tsamp']
period = np.float128(pulsar_period) / tsamp
window = int(period)




#offset = period-window
#ind = []
#for i in range(1,data.shape[0]/window):
  #offset += offset
  #if offset > 1: 
    #ind.append(window*i)
    #offset -= 0.9999999999


ind = np.arange(data.shape[0]/window)*period
ind = ind.astype(np.int)%window
ind = np.where(np.abs(np.diff(ind))>0)[0]*window+window
ind += np.arange(ind.size)




data = np.delete(data,ind,axis=0)
data = data[:data.shape[0]/window*window]
data = np.mean(data.reshape(-1, window, nchans),axis=0)
ind = np.arange(0,2592,16)
data = np.delete(data,ind,axis=1)


#Plot 

def color_range(data):
  #Define the color range
  clean = data[data>0]
  min_element = clean.size/20
  max_element = clean.size*9/10
  vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
  vmax = np.partition(clean, max_element, axis=None)[max_element]
  return vmin,vmax

vmin,vmax = color_range(data)

plt.figure(figsize=(20,10))
plt.xlabel('Time (s)')
plt.ylabel('Frequency (MHz)')
plt.imshow(data.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest',extent=[t0,t0+total_duration,119,151],vmin=vmin,vmax=vmax)
plt.colorbar()
x_lines = np.linspace( t0+res, t0+res*(n_bins-1), n_bins-1 )
y_lines = np.zeros(x_lines.shape[0])
plt.plot([x_lines,x_lines],[y_lines+118,y_lines+151],'b--')
plt.axis([t0,t0+total_duration,119,151])

#plot the DM lines
DM=300
time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
plt.plot(time,freq,'b-')
plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')

