#############################
#
# Parameters
#
# Written by Daniele Michilli
#
#############################


DM_MIN = 5.  #pc/cm3
SIGMA_MIN = 10.
#DURATION_MAX = 0.049152  #s


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz

ERR_FLOAT = 0.001
STEPS_GROUP = 3
DURAT_GROUP = 1

SIGMA_TOLL = 4
SIGMA_TOLL_IB = 2

RFI_percent = 3
PULS_LENGHT = 15

FILTERS = [
  [3.44e-2,3.44e-2,3.44e-2],    #duration
  [0.9,0.773,0.675],            #dDM/(N_events-1)/step
  [1.,1.,1.],                   #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.03,1.05,1.08],             #puls.Sigma/Sigma_min
  [6.33e-3,7.52e-3,3.41e-4],    #puls.Sigma/Sigma_min**4
  [1.45,1.3,10.2],              #abs(Sigma_DM_max-Sigma_DM_min)
  [3.04,9.06,2.83],             #f (>)
  [-41.,-36.3,-48.8],           #f (<)
  [3.61e-5,8.24e-5,1.28e-3]]    #gb.Time.apply(np.var)

FILTERS_BEST = [
  [0.667,0.636,0.537],          #dDM/(N_events-1)/step
  [0.667,0.636,0.714],          #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.06,1.09,1.19],             #puls.Sigma/Sigma_min
  [8.06e-3,8.36e-3,9.26e-3],    #puls.Sigma/Sigma_min**4
  [0.54,0.41,3.39],             #abs(Sigma_DM_max-Sigma_DM_min)
  [-12.2,-14.2,-14.6],          #f (>)
  [-25.7,-26.1,-27.1],          #f (<)
  [3.09e-6,7.28e-6,6.34e-5]]    #gb.Time.apply(np.var)

FILTERS_INC = [
  [7.37e-02,7.37e-02],  #duration
  [8.33e-01,5.67e-01],  #dDM/(N_events-1)/step
  [1.,7.60e-01],        #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.02,1.10],          #puls.Sigma/Sigma_min
  [22.8,2.16],          #abs(Sigma_DM_max-Sigma_DM_min)
  [-41.0,-69.1],        #puls.dDM/puls.dTime
  [8.52e-06,5.61e-03]]  #gb.Time.apply(np.var)

FILTERS_BEST_INC = [
  [6.43e-01,5.15e-01],  #dDM/(N_events-1)/step
  [1.,3.85e-01],        #abs(puls.DM-puls.DM_c)/puls.dDM
  [1.05,1.22],          #puls.Sigma/Sigma_min
  [0.34,0.29],          #abs(Sigma_DM_max-Sigma_DM_min)
  [-20.5,-26.1],        #puls.dDM/puls.dTime
  [3.17e-06,1.18e-03]]  #gb.Time.apply(np.var)
