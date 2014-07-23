cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

def group(double[::1] DM not None,
          double[::1] Sigma not None,
          double[::1] Time not None,
          double[::1] Down_Time not None,
          long[::1] Pulse not None):
  
  
  cdef:
    unsigned int i, j, j_min, empty, SNR_max
    unsigned int n_steps = 2
    unsigned int j_max = 0
    unsigned int code = 0
    unsigned int dim = len(DM)
    double step_min, step_max, dDM
    double DM_new=-1.
    double float_err=0.0001


  
  for i in range(0,dim):
    
    if Pulse[i]==-1: continue
    
    if Pulse[i]==0: 
      code += 1
      Pulse[i] = code
  
    if DM[i] != DM_new:
      
      j_min = 0
      
      DM_new = DM[i]
      
      if DM_new < 40.: 
        
        step_min = 0.01-float_err
        
        step_max = n_steps*0.01+float_err
      
      elif DM_new < 140.: 
        
        step_min = 0.05-float_err
        
        step_max = n_steps*0.05+float_err
      
      else: 
        
        step_min = 0.1-float_err
        
        step_min = n_steps*0.1+float_err
        
      for j in range(i+1,dim):
      
        dDM = DM[j] - DM_new

        if (dDM>step_min) and (dDM<step_max): 
          
          if j_min == 0: j_min = j
          j_max = j
          
    
    empty = 0
          
    for j in range(j_min,j_max):
      
      if abs(Time[i]-Time[j]) < Down_Time[i]+Down_Time[j]:
        
        if empty == 0:
          
          Pulse[j] = code
          SNR_max = j
          empty = 1
          
        else:
          
          if Sigma[j] > Sigma[SNR_max]:
            
            Pulse[SNR_max] = -1
            SNR_max = j
            Pulse[j] = code
            
          else:
            
            Pulse[j] = -1
  
  return