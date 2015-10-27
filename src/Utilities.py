import numpy as np
import pandas as pd
import os
import math
import logging
import sys

try:
  import filterbank
  import sigproc
  import pyfits
except ImportError: pass

from Parameters import *


def read_filterbank(filename,DM,bin_start):
  if (not 'filterbank' in sys.modules) or (not 'sigproc' in sys.modules):
    logging.warning("Additional modules missing")
    return
  
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
  #ind = np.arange(0,2592,16)
  #spectrum = np.delete(spectrum,ind,axis=1)
  
  return spectrum
 


def read_fits(fits,DM,bin_start,duration,offset,RFI_reduct=False):
  bin_start = np.int(bin_start)
  duration = np.int(duration/RES+1)
  
  if isinstance(fits,str):
    try: fits = pyfits.open(fits,memmap=True)
    except NameError: 
      logging.warning("Additional modules missing")
      return    
    header = fits['SUBINT'].header
  else: header = fits['SUBINT'].header
  
  N_channels = header['NCHAN']
  N_spectra = header['NSBLK']

  bin_start -= offset
  bin_end = bin_start + DM2delay(DM) + duration + 2*offset
  
  subint_start = bin_start/N_spectra
  subint_end = bin_end/N_spectra+1
  subint = fits['SUBINT'].data[subint_start:subint_end]['DATA']
  
  fits.close()

  subint = subint.reshape((subint_end-subint_start)*N_spectra,N_channels)
  
  #ATTENZIONE! Salvare .mask file nella pipeline!!
  try:
    #Zap the channels from the mask file
    zap_chans, zap_ints, chans_per_int = read_mask(filename,subint_start,subint_end)
    med_value = np.median(subint)
    subint[:,zap_chans] = med_value
    zap_ints = zap_ints[(zap_ints>=bin_start)&(zap_ints<=bin_end)]
    subint[zap_ints,:] = med_value

    for i,chans in enumerate(chans_per_int):
      chunk = subint[i*N_spectra:(i+1)*N_spectra]
      chunk[:,chans] = med_value
      if RFI_reduct: np.clip(chunk - np.median(chunk,axis=0) + 128, 0, 255, out=chunk)
    
  except IOError:
    if RFI_reduct:
      for i in range(subint_end-subint_start):
        chunk = subint[i*N_spectra:(i+1)*N_spectra]
        np.clip(chunk - np.median(chunk,axis=0) + 128, 0, 255, out=chunk)
    
  subint = subint[bin_start%N_spectra:bin_start%N_spectra+bin_end]
  #ind = np.arange(0,2592,16)
  #subint = np.delete(subint,ind,axis=1)
  
  return subint



def read_mask(filename,subint_start,subint_end):
  f = open(filename,'r')
  timesigma, freqsigma, mjd, dtint, lofreq, dfreq = np.fromfile(f,count=6,dtype=np.double)
  numchan, numint, ptsperint = np.fromfile(f,count=3,dtype=np.intc)
  
  num_zap_chans, = np.fromfile(f,count=1,dtype=np.intc)
  zap_chans = np.fromfile(f,count=num_zap_chans,dtype=np.intc)

  num_zap_ints, = np.fromfile(f,count=1,dtype=np.intc)
  zap_ints = np.fromfile(f,count=num_zap_ints,dtype=np.intc)
  
  num_chans_per_int = np.fromfile(f,count=numint,dtype=np.intc)
  num_chans_per_int[num_chans_per_int==numchan] = 0
  num_chans_cum = np.cumsum(num_chans_per_int)
  chans = np.fromfile(f,dtype=np.intc)
  f.close()
  
  chans_per_int = []
  for i in range(subint_start,subint_end+1):
    if i==0:
      chans_per_int.append(chans[:num_chans_cum[i]])
    else:
      chans_per_int.append(chans[num_chans_cum[i-1]:num_chans_cum[i]])
    
  return zap_chans, zap_ints, chans_per_int
  


def clean(data):

  data_wso = data.astype(np.int)

  masked = np.ma.array(data_wso)
  masked[:,np.arange(0,2592,16)] = np.ma.masked


  #Set the mean of each spectra to zero
  #med_col = np.mean(masked,axis=1)
  med_col = np.min((\
   np.mean(masked[:,:2592/5],axis=1),\
   np.mean(masked[:,-2592/5:],axis=1)),\
   axis=0)

  data_wso = data_wso - med_col[:,np.newaxis]

  #Shift the mean to 128
  data_wso = data_wso + 128

  #Set the right proprieties to the data
  data_wso[:,np.arange(0,2592,16)] = 0
  dtype_min = np.iinfo(data.dtype).min
  dtype_max = np.iinfo(data.dtype).max
  np.clip(data_wso, dtype_min, dtype_max, out=data_wso)
  data_wso = np.around(data_wso)
  data = data_wso.astype(data.dtype)

  return data

 
 

def read_header(filename):
  name, ext = os.path.splitext(filename)
  if ext == '.fits': 
    if not 'pyfits' in sys.modules:
      logging.warning("Additional modules missing")
      return    
    fits = pyfits.open(filename,memmap=True)
    header = fits['SUBINT'].header + fits['PRIMARY'].header
    fits.close()
    return header
  elif ext == '.fil':
    if not 'filterbank' in sys.modules:
      logging.warning("Additional modules missing")
      return
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

#Probability that k elements of an ensamble n will group together in a space with n dimensions
#Formula approssimata, ricalcolare quella esatta
def p(n,k,dim = 5450.):
  return c(n,k)/dim**(k-1)



def d(n,k):                           
  return math.factorial(n)/math.factorial(n-k)
#def p(n,k,dim):
  #return 1-float(d(dim,(n-k+2)))*dim**(k-2)/dim**n
#float(Utilities.c(n,k)*dim**(n-k)-n+1)*dim/dim**n




#test of the probability formula
def test_p(n,k,tot,dim = 5450.):
  dim = float(dim)
  fav = 0
  for i in range(tot):
    test = pd.Series(np.random.randint(0,dim,n))
    test = test.groupby(by=test).size().max()
    fav += test>=k 
  return np.float(fav)/tot

def color_range(data):
  #Define the color range
  clean = data[data>0]
  min_element = clean.size/20
  max_element = clean.size*9/10
  vmin = np.partition(clean, min_element, axis=None)[min_element]   #invece di subint[subint>0] possibile subint[:-(num_rows/down_fact)]
  vmax = np.partition(clean, max_element, axis=None)[max_element]
  return vmin,vmax


def rrat_period(times, numperiods=20000):
  #Modified version of PRESTO
  ts = np.asarray(sorted(times))
  ps = (ts[1]-ts[0])/np.arange(1, numperiods+1)
  dts = np.diff(ts)
  xs = dts / ps[:,np.newaxis]
  metric = np.sum(np.fabs((xs - xs.round())), axis=1)
  pnum = metric.argmin()
  numrots = xs.round()[pnum].sum()
  p = (ts[-1] - ts[0]) / numrots
  rotations = dts / p
  diff_max = np.max(np.abs(np.round(rotations) - rotations))
  return p, diff_max
  
  
def DynamicSpectrum(pulses,idL,sap,beam,store):    #Creare una copia di pulses quando si chiama la funzione!
  plt.clf()
  
  if beam==12: stokes = 'incoherentstokes'
  else: stokes = 'stokes'
  filename = '{folder}/{idL}_red/{stokes}/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}.fits'.format(folder=Paths.RAW_FOLDER,idL=idL,stokes=stokes,sap=sap,beam=beam)
  if not os.path.isfile(filename): return
  
  
  pulses.Sample[pulses.DM>141.71] *= 4
  pulses.Sample[(pulses.DM>40.47)&(pulses.DM<=141.71)] *= 2
  
  #controllare che questo vada dopo downsampling correction!
  header = Utilities.read_header(filename)
  MJD = header['STT_IMJD'] + header['STT_SMJD'] / 86400.
  v = presto.get_baryv(header['RA'],header['DEC'],MJD,1800.,obs='LF')
  bin_idx = np.int(np.round(1./v))
  pulses.Sample += pulses.Sample/bin_idx
  
  
  
  #fig = plt.figure()
  
  freq = np.arange(151,118,-1,dtype=np.float)
  offset = 300
  
  #for i,(idx,puls) in enumerate(pulses.iterrows()):
    #spectrum, bin_start, bin_end = read_fits(filename,puls.DM,puls.Sample,offset)
    ##if RFIexcision.Multimoment() > FILTERS['Multimoment']: 
      ##pulses.Pulse += 1
      ##continue
    #extent = [bin_start*RES,bin_end*RES,119,151]
    #ax = plt.subplot2grid((2,5),(i/5,i%5))
    #ax.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent)
    #time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2)
    #ax.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
    ##ax.plot((time[0],time[0]+puls.Sample*RES),(freq[0],freq[0]),'r',(time[-1],time[-1]+puls.Sample*RES),(freq[-1],freq[-1]),'r')
    #ax.axis(extent)
    #ax.set_title('Sigma = {0:.1f}, Rank = {1}'.format(puls.Sigma,i))



  puls = pulses
  spectrum, bin_start, bin_end = read_fits(filename,puls.DM,puls.Sample,offset)
  extent = [bin_start*RES,bin_end*RES,119,151]
  
  vmin,vmax = color_range(spectrum)
  
  plt.imshow(spectrum.T,cmap='Greys',origin="lower",aspect='auto',interpolation='nearest',extent=extent,vmin=vmin,vmax=vmax)     #provare anche pcolormesh
  time = 4149 * puls.DM * (np.power(freq,-2) - 151.**-2) + bin_start*RES
  plt.plot(time-offset*RES,freq,'r',time+offset*RES,freq,'r')
  plt.axis(extent)
  
  
  
  
  
  #allineamento perfetto
  #problema di contrasto: provare a diminuire range di colori, azzerare minori di soglia, mediare in duration del pulse, etc
  #prova di de-dispersione (forse inutile):
  freq = np.linspace(F_MIN,F_MAX,2430)
  time = (4149 * DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
  for i in range(spectrum.shape[1]):
    spectrum[:,i] = np.roll(spectrum[:,i], time[i])

  
  
  
  
  # Set common labels
  #fig.text(0.5, 0.05, 'Time (s)', ha='center', va='center', fontsize=8)
  #fig.text(0.08, 0.5, 'DM (pc/cm3)', ha='left', va='center', rotation='vertical', fontsize=8)
  #fig.text(0.5, 0.95, str(idL), ha='center', va='center', fontsize=12)
  
  #plt.savefig(store,format='png',bbox_inches='tight',dpi=200)
  
  return 





def dedisp(spectrum,DM,duration):
  freq = np.linspace(F_MIN,F_MAX,spectrum.shape[1])
  time = (4149 * DM * (F_MAX**-2 - np.power(freq,-2)) / RES).round().astype(np.int)
  for i in range(spectrum.shape[1]):
    spectrum[:,i] = np.roll(spectrum[:,i], time[i])
  duration = np.int(np.round(duration/RES))
  
  
