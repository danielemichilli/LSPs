#############################
#
# Cython Functions
#
# Functions written in cython
# are grouped here.
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
    float durat = DURAT_GROUP
    unsigned int code = 0
    unsigned int dim = len(DM)
    float step, step_min, step_max, dDM, DM_min
    float DM_new = -1.
    float float_err = 0.001


  # Assign a code to each event.
  # Events must have been ordered in DM and then in Time
  for i in range(0,dim):
    
    #Remove close events at the same DM
    j = i+1
    if DM[i] == DM[j]:
    
      if abs(Time[i]-Time[j]) < durat:
        
        if j < dim : 
        
          if Sigma[i] < Sigma[j] : Pulse[i] = -1
          else : Pulse[j] = -1
  
    if Pulse[i]==-1: continue
    
    # Set a code to the events that aren't in a pulse
    if Pulse[i]==0: 
      code += 1
      Pulse[i] = code
      
    # Defines for a certain DM a range of events that can be grouped
    if DM[i] != DM_new:
      
      j_min = 0
      j_max = dim
      
      DM_new = DM[i]
      
      if DM_new < 40.49: step = 0.01
        
      elif DM_new < 141.69: step = 0.05
        
      else: step = 0.1
        
      step_min = step - float_err
      
      step_max = step * (n_steps + 1) + float_err
      
      
      #find the minimum and maximum event in the range
      for j in range(i+1,dim):
        
        dDM = DM[j] - DM_new

        if dDM > step_max:
          
          j_max = j
          
          break

        if dDM > step_min: 
          
          if j_min == 0: j_min = j
          
    empty = 0
    
    if j_min > 0:

      # Gives a code to the next event in the pulse
      for j in range(j_min,j_max):

        if abs(Time[i]-Time[j]) < durat:   #MAYBE add a condition on SNR (attention: dSNR depends on SNR!) 
          
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
    unsigned int i, j, RFI
    unsigned int dim_l = len(DM_c_l)
    unsigned int dim_r = len(DM_c_r)
    unsigned int j_min = 0
    int TollSigma
    int rfi_limit = RFI_percent
    float DTime, Time, DM, DDM, sign
    float Duration_Max = FILTERS[0][2]
  
  # Uses different tollerances on sigma for coherent and incoherent beams
  if CB==int(1): TollSigma = SIGMA_TOLL
  else: TollSigma = SIGMA_TOLL_IB

  # Compare each pulse of the first group with each of the second
  # Assign different codes under certain conditions 
  for i in range(0,dim_l):
    
    RFI = 0
    if Pulse_l[i] >= rfi_limit: RFI = 1 
    
    j_flag = j_min
    
    for j in range(j_min, dim_r):
      
      if Pulse_r[j] >= rfi_limit: 
        
        if RFI > 0:
        
          continue
      
      Time = abs(Time_c_l[i]-Time_c_r[j])
      DTime = dTime_l[i]+dTime_r[j]
      
      sign = Time_c_l[i]-Time_c_r[j]
      
      if Time < 4.*DTime :
      
        if j_flag == j_min: j_min = j
        
        DM = abs(DM_c_l[i]-DM_c_r[j])
        DDM = dDM_l[i]+dDM_r[j]
        
        if DM < 2.*DDM :
          
          if abs(Sigma_l[i]-Sigma_r[j]) < TollSigma:
                        
            Pulse_l[i] += 1
            Pulse_r[j] += 1
          
      elif sign < 0.: 
        
        break
      
  return
  
  
  
  
  
  
  
  
#-------------------------
# Compares repeated pulses
#-------------------------
def Compare_candidates(float[::1] DM not None,
            float[::1] Sigma not None,
            float[::1] Time not None,            
            short[::1] N_pulses not None,
            long[::1] idx not None,            
            long[::1] cand not None):

  cdef:
    unsigned int i, j
    unsigned int dim = len(DM)

  # Compare each candidate
  for i in range(0,dim):
    
    for j in range(dim-1,i,-1):
      
      if abs(DM[j]-DM[i]) < 1.:
      
        if abs(Time[j]-Time[i]) < 1.:
        
          cand[i] = idx[j]
          
          break

      
  return