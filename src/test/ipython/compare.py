%%cython
import pandas as pd
import numpy as np

cimport numpy as np
from cpython cimport bool



SIGMA_TOLL = 4


idL = 'L98709'
store = pd.HDFStore('/home/michilli/Documents/L98709/SinlgePulses.hdf5','r')
p0 = store[idL+'_pulses']
store.close()

idL = 'L196710'
store = pd.HDFStore('/home/michilli/Documents/L196710/SinlgePulses.hdf5','r')
p1 = store[idL+'_pulses']
store.close()




#data.sort(['SAP','BEAM','DM'],inplace=True)

puls0 = p0.ix[:,['DM','dDM','Time','Duration','Sigma']].astype(np.float32)
puls1 = p1.ix[:,['DM','dDM','Time','Duration','Sigma']].astype(np.float64)


cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

def Compare_CB(float[::1] DM_l not None,
               float[::1] dDM_l not None,
               float[::1] Time_l not None,
               float[::1] Duration_l not None,
               float[::1] Sigma_l not None,
               int[::1] Pulse_l not None,
               double[::1] DM_r not None,
               double[::1] dDM_r not None,
               double[::1] Time_r not None,
               double[::1] Duration_r not None,
               double[::1] Sigma_r not None,
               long[::1] Pulse_r not None):

  cdef:
    unsigned int i, j
    unsigned int dim_l = len(DM_l)
    unsigned int dim_r = len(DM_r)
    unsigned int TollSigma = SIGMA_TOLL


  for i in range(0,dim_l):
        
    for j in range(0,dim_r):
      
      if abs(Time_l[i]-Time_r[j]) < Duration_l[i]+Duration_r[j]:  #Duration cattivo parametro: RFI hanno Sigma_max casuale, meglio dTime
        
        if abs(DM_l[i]-DM_r[j]) < dDM_l[i]+dDM_r[j]:  #Condizione forte: pulses eliminati anche se c'e' overlapping minimo
          
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
            
            Pulse_l[i] = 0
            Pulse_r[j] = 0

  return

p0['Pulse']=1
p0['Pulse']=p0.Pulse.astype(np.intc)
p1['Pulse']=1

Compare_CB(puls0.DM.values,puls0.dDM.values,puls0.Time.values,puls0.Duration.values,puls0.Sigma.values,p0.Pulse.values,\
           puls1.DM.values,puls1.dDM.values,puls1.Time.values,puls1.Duration.values,puls1.Sigma.values,p1.Pulse.values)
  
  
#store = pd.HDFStore('SinlgePulses.hdf5','w')
#store[idL] = data
#store.close()

