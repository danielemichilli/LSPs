import numpy as np
import os
import math
import logging

try:
  import filterbank
  import sigproc
  import pyfits
except ImportError: pass

from Parameters import *


def read_filterbank(filename,DM,bin_start):
  try:  #provare metodi migliori
    id(filterbank)
    id(sigproc)
  except NameError: 
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
  ind = np.arange(0,2592,16)
  spectrum = np.delete(spectrum,ind,axis=1)
  
  return spectrum
 


def read_fits(filename,DM,bin_start,offset):
  try:  #provare metodi migliori
    id(pyfits)
  except NameError: 
    logging.warning("Additional modules missing")
    return

  bin_start = np.int(bin_start)
  fits = pyfits.open(filename,memmap=True)
  
  header = fits['SUBINT'].header
  
  N_channels = header['NCHAN']
  N_spectra = header['NSBLK']

  bin_start -= offset
  bin_end = bin_start + DM2delay(DM) + 2*offset
  
  subint_start = bin_start/N_spectra
  subint_end = bin_end/N_spectra+1
  subint = fits['SUBINT'].data[subint_start:subint_end]['DATA']
  
  fits.close()

  subint = subint.reshape((subint_end-subint_start)*N_spectra,N_channels)
  
  subint = subint[bin_start%N_spectra:bin_start%N_spectra+bin_end]
  ind = np.arange(0,2592,16)
  subint = np.delete(subint,ind,axis=1)
  
  return subint
 
 

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

#Probability that k elements of an ensamble n will group together in a space with n dimensions
def p(n,k):
  dim = 5450.
  return c(n,k)/dim**(k-1)

#test of the probability formula
def test_p(n,k):
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
    #spectrum, bin_start, bin_end = Utilities.read_fits(filename,puls.DM,puls.Sample,offset)
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
  spectrum, bin_start, bin_end = Utilities.read_fits(filename,puls.DM,puls.Sample,offset)
  extent = [bin_start*RES,bin_end*RES,119,151]
  
  vmin,vmax = Utilities.color_range(spectrum)
  
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
  
  
