#!/usr/bin/env python

'''

Downsampler

Written by Daniele Michilli

'''


import numpy as np
import os
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits


def dynspect_plot(idx,idL,filename,raw_file):
  
  #Constants declaration
  res = 36 #s  #bin resolution
  n_bins = 5  #number of bins to plot
  num_of_spectra = 1000  #number of spectra to actually plot

  #Raw data file
  fits = pyfits.open(raw_file,memmap=True)

  #Parameters for the spectrum
  subint_duration = fits['SUBINT'].header['TBIN']*fits['SUBINT'].header['NSBLK']
  total_duration = res*n_bins #s 
  frame_duration = np.int(total_duration/subint_duration)
 
  #Downsample the data
  N_channels = fits['SUBINT'].header['NCHAN']
  N_spectra = frame_duration*fits['SUBINT'].header['NSBLK']
  down_fact = N_spectra / num_of_spectra
  while N_spectra % down_fact != 0: down_fact -= 1 #find the closest integer divisor to average

  #Prepare the DM lines to plot
  freq = np.arange(151,117,-1,dtype=np.float)
  DM=300
  time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
  DM=500
  time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1
  DM=700
  time = 4149 * DM * (np.power(freq,-2) - 151.**-2) + t0 + 1


  for ind0 in idx:
    #Set the start of the spectrum
    t0 = (ind0-1) * res
    subint_index = np.int(t0/subint_duration)

    #Load the data
    subint = fits['SUBINT'].data[subint_index:subint_index+frame_duration]['DATA']

    #Average and clean the spectrum
    subint = subint.reshape(N_spectra/down_fact,down_fact,N_channels).mean(axis=1)
    subint = subint[:,~np.all(subint == 0, axis=0)]

    #Define the color range
    min_element = subint.size/20
    max_element = subint.size*9/10
    vmin = np.partition(subint, min_element, axis=None)[min_element]
    vmax = np.partition(subint, max_element, axis=None)[max_element]
    
    #Plot the spectrum
    plt.figure(figsize=(20,10))
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (MHz)')
    plt.imshow(subint.T,cmap=mpl.cm.hot_r,origin="lower",aspect='auto',interpolation='nearest',extent=[t0,t0+total_duration,119,151],vmin=vmin,vmax=vmax)
    plt.colorbar()
    x_lines = np.linspace( t0+res, t0+res*(n_bins-1), n_bins-1 )
    y_lines = np.zeros(x_lines.shape[0])
    plt.plot([x_lines,x_lines],[y_lines+118,y_lines+151],'b--')
    plt.axis([t0,t0+total_duration,119,151])

    #plot the DM lines
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')
    plt.plot(time,freq,'b-')
    plt.annotate(str(DM),(time[freq==120],120),color='b',horizontalalignment='left',fontweight='bold')

    plt.savefig('{}_{}_DynSpect_{}.pnx'.format(idL,filename,ind0),format='png',bbox_inches='tight',dpi=150)

  fits.close()




print 'Starting to downsample'

#Define the parser
parser = argparse.ArgumentParser()
parser.add_argument('idL', nargs=1)
parser.add_argument('filenm', nargs=1)
parser.add_argument('raw_file', nargs=1)
args = parser.parse_args()
filenm = args.filenm[0]
idL = args.idL[0]
raw_file = args.raw_file[0]

#Load the timeserie
timeserie = np.fromfile('{}_{}.dat'.format(idL,filenm),dtype=np.float32)

#Set all the timeseries to start at the same time   ---   CHECK THIS!!
if filenm != 'DM0.00':
  #Check the starting time of the current DM timeseries
  inf_file = open('{}_{}.inf'.format(idL,filenm))
  for line in inf_file:
    if line.startswith(" Epoch of observation"):
      mjd = line.split('=  ')[-1]
    #if line.startswith(" Width of each time series bin"):
      #res = line.split('=  ')[-1]
  inf_file.close()

  #Check the starting time of DM0 timeseries
  inf_file = open('{}_DM0.00.inf'.format(idL))
  for line in inf_file:
    if line.startswith(" Epoch of observation"):
      mjd_DM0 = line.split('=  ')[-1]
  inf_file.close()
  
  #Calculate the time delay
  dt = (float(mjd_DM0) - float(mjd)) *86400 #s
  res = 0.000491519982460886*4
  didx = int(dt/res)
  
  #Apply the time delay
  timeserie = np.roll(timeserie,didx)

#Downsample to 36s time resloution
if timeserie.size == 7392000:  down_fact = 73920
elif timeserie.size == 3696000: down_fact = 36960
elif timeserie.size == 1848000: down_fact = 18480
else:
  print "Error: length of the timeserie unknown!"
  exit()

downsampled = np.mean(timeserie.reshape(100,down_fact),axis=1)  #The first value is the length after the downsample, the second the downsampling factor
downsampled.tofile('{}_{}_down_30s.ds'.format(idL,filenm))

#Plot bins that are some factor above the median and that have higher signal compared to DM0
if filenm != 'DM0.00':
  med = np.median(downsampled)
  DM0 = np.fromfile('{}_DM0.00_down_30s.ds'.format(idL),dtype=np.float32)
  idx = np.where( ( downsampled > 1.005 * med )&( downsampled > DM0 )&( downsampled > np.roll(DM0,-1) )&( downsampled > np.roll(DM0,-2) ))[0]
  idx.tofile('{}_{}_down_30s(ind_cand).dx'.format(idL,filenm))
  dynspect_plot(idx,idL,filenm,raw_file)
  
print 'Finished to downsample' 

