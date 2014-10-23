%%cython
import pandas as pd
import numpy as np

cimport numpy as np
from cpython cimport bool


idL = 'L196710'
sap = 1
beam = 71
'''
idL = 'L98709'
sap = 0
beam = 54
'''


store = pd.HDFStore('_SinlgePulses.hdf5','r')
data = store[idL]
store.close()


#data = data[(data['DM']>43) & (data['DM']<44) & (data['Time']>844) & (data['Time']<848)]



data.sort(['SAP','BEAM','DM'],inplace=True)


data['Pulse']=0

data_copy = data.ix[:,['DM','Sigma','Time','Down_Time']].astype(np.float64)



cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

def group(double[::1] DM not None,
          double[::1] Sigma not None,
          double[::1] Time not None,
          double[::1] Down_Time not None,
          long[::1] Pulse not None):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
    unsigned int n_steps = 4
    unsigned int code = 0
    unsigned int dim = len(DM)
    double step, step_min, step_max, dDM, DM_min
    double DM_new = -1.
    double float_err = 0.0001


  
  for i in range(0,dim):
    
    if Pulse[i]==-1: continue
    
    if Pulse[i]==0: 
      code += 1
      Pulse[i] = code
      
  
    if DM[i] != DM_new:
      
      j_min = 0
      j_max = dim
      
      DM_new = DM[i]
      
      if DM_new < 40.: step = 0.01
        
      elif DM_new < 140.: step = 0.05
        
      else: step = 0.1
        
      step_min = step - float_err
      
      step_max = step * n_steps + float_err
        
        
      for j in range(i+1,dim):
        
        dDM = DM[j] - DM_new

        if dDM > step_max:
          
          j_max = j
          
          break

        if dDM > step_min: 
          
          if j_min == 0: j_min = j
        

    empty = 0
    
    if j_min > 0:
      
      for j in range(j_min,j_max):
        
        if abs(Time[i]-Time[j]) < Down_Time[i]+Down_Time[j]:
          
          if empty == 0:
            
            Pulse[j] = Pulse[i]
            SNR_max = j
            empty = 1
            DM_min = DM[j]
            
          else:
            
            if DM[j] > DM_min: break
            
            if Sigma[j] > Sigma[SNR_max]:
              
              Pulse[SNR_max] = -1
              SNR_max = j
              Pulse[j] = Pulse[i]
              
            else:
              
              Pulse[j] = -1
    
  return

group(data_copy.DM.values,data_copy.Sigma.values,data_copy.Time.values,data_copy.Down_Time.values,data.Pulse.values)

#data = data[data.Pulse>0]

store = pd.HDFStore('SinlgePulses.hdf5','w')
store[idL] = data
store.close()

