cimport cython

#---------------------------------
# Gives a pulse-code to each event
#---------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Get_Group(double[::1] DM not None,
          double[::1] Sigma not None,
          double[::1] Time not None,
          double[::1] Duration not None,
          int[::1] Pulse not None,
          double step,
          double durat,
          double dDM_max):
  
  cdef:
    unsigned int i, j, k, j_min, j_max, empty, SNR_max
    unsigned int n_steps = int(dDM_max / step)
    unsigned int code = 0
    unsigned int dim = len(DM)
    double step_min, step_max, dDM, DM_min
    double DM_new = -1.
    double float_err = 0.000001


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

        if abs(Time[i]-Time[j]) < durat:
          
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
   

   
   
#-------------------------------------
# Gives a candidate-code to each pulse
#-------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def Get_Candidate(
          long[::1] Beam not None,
          double[::1] Sigma not None,
          double[::1] Time not None,
          int[::1] Pulse not None,
          double[::1] Duration not None):
  
  cdef:
    unsigned int dim = len(Time)
    unsigned int code = 0
    
  for i in range(0,dim):
  
    if Pulse[i] == -1: continue

    if Pulse[i] == -2:
      
      Pulse[i] = code
      
      code += 1
    
    else:
      
      for j in range(i+1,dim):
        
        if Time[j] - Time[i] < Duration * 2.: 
          
          if Beam[i] != Beam[j]: Pulse[j] = code
            
          elif Sigma[i] > Sigma[j]: Pulse[j] = -1
            
          else: Pulse[i] = -1
                  
        else: break
         
  return