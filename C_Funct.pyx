#############################
#
# Cython Functions
#
# Written by Daniele Michilli
#
#############################

#Installation: sudo python setup.py build_ext --inplace

cimport cython
from Parameters import *

@cython.boundscheck(False)
@cython.wraparound(False)

#---------------------------------
# Gives a pulse-code to each event
#---------------------------------
def Get_Group(float[::1] DM not None,
          float[::1] Sigma not None,
          float[::1] Time not None,
          float[::1] Duration not None,
          int[::1] Pulse not None):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
    unsigned int n_steps = STEPS_GROUP
    unsigned int durat = DURAT_GROUP
    unsigned int code = 0
    unsigned int dim = len(DM)
    float step, step_min, step_max, dDM, DM_min
    float DM_new = -1.
    float float_err = 0.0001


  
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
        
        if abs(Time[i]-Time[j]) < durat*(Duration[i]+Duration[j]):
          
          if Pulse[j] == -1: continue
          
          if Pulse[j] > 0: 
            
            Pulse[j] = -1
            continue
          
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
  
  
#-------------------------
# Compares different beams
#-------------------------
def Compare(float[::1] DM_c_l not None,
            float[::1] dDM_l not None,
            float[::1] Time_c_l not None,
            float[::1] dTime_l not None,
            float[::1] Sigma_l not None,
            signed char[::1] Pulse_l not None,
            float[::1] DM_c_r not None,
            float[::1] dDM_r not None,
            float[::1] Time_c_r not None,
            float[::1] dTime_r not None,
            float[::1] Sigma_r not None,
            signed char[::1] Pulse_r not None,
            int CB):

  cdef:
    unsigned int i, j
    unsigned int dim_l = len(DM_c_l)
    unsigned int dim_r = len(DM_c_r)
    int TollSigma
    float DTime,Time,DM,DDM
    
  if CB==int(1): TollSigma = SIGMA_TOLL
  else: TollSigma = SIGMA_TOLL_IB

  for i in range(0,dim_l):
        
    for j in range(0,dim_r):
      
      Time = abs(Time_c_l[i]-Time_c_r[j])
      DTime = dTime_l[i]+dTime_r[j]
      
      if Time < DTime :
      
        DM = abs(DM_c_l[i]-DM_c_r[j])
        DDM = dDM_l[i]+dDM_r[j]
        
        if DM < DDM :
          
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
            
            Pulse_l[i] += 2
            Pulse_r[j] += 2
          
        elif DM < 10.*DDM :
        
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
          
            Pulse_l[i] += 1
            Pulse_r[j] += 1
      
      elif Time < 6. * DTime :
      
        if abs(DM_c_l[i]-DM_c_r[j]) < 10.*(dDM_l[i]+dDM_r[j]) :
        
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
          
            Pulse_l[i] += 1
            Pulse_r[j] += 1
      
  return