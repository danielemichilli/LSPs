import numpy as np
import pandas as pd

import C_Funct
import RFIexcision
from Parameters import *


def Generator(events):
  #-------------------------------
  # Create a table with the pulses
  #-------------------------------
    
  gb = events.groupby('Pulse',sort=False)
  pulses = events.loc[events.index.isin(gb.Sigma.idxmax())]
  pulses.index = pulses.Pulse
  pulses.index.name = None
  pulses = pulses.loc[:,['SAP','BEAM','DM','Sigma','Time','Duration']]
  pulses['Pulse'] = 0
  pulses.Pulse = pulses.Pulse.astype(np.int8)
  pulses['Candidate'] = -1
  pulses.Candidate = pulses.Candidate.astype(np.int8)
  pulses['dDM'] = (gb.DM.max() - gb.DM.min()) / 2.
  pulses.dDM=pulses.dDM.astype(np.float32)
  pulses['dTime'] = (gb.Time.max() - gb.Time.min()) / 2.
  pulses.dTime=pulses.dTime.astype(np.float32)
  
  #pulses['dSample'] = (gb.Sample.max() - gb.Sample.min()) / 2.
  #pulses.dSample = pulses.dSample.astype(np.float32)  
  
  pulses['DM_c'] = (gb.DM.max() + gb.DM.min()) / 2.
  pulses.DM_c=pulses.DM_c.astype(np.float32)
  pulses['Time_c'] = (gb.Time.max() + gb.Time.min()) / 2.
  pulses.Time_c=pulses.Time_c.astype(np.float32)
  pulses['N_events'] = gb.DM.count()
  pulses.N_events = pulses.N_events.astype(np.int16)

  #####   MODIFICARE   #########
  
  pulses = pulses[pulses.N_events>4]
  
  # Reduce the RFI and corrects for the time misalignment
  
  #data_idxmax = data.loc[gb.DM.idxmax()]
  #data_idxmax.index = data_idxmax.Pulse
  #data_idxmin = data.loc[gb.DM.idxmin()]
  #data_idxmin.index = data_idxmin.Pulse  
  #Sigma_DM_max = data_idxmax.Sigma
  #Sigma_DM_min = data_idxmin.Sigma
  #Time_DM_max = data_idxmax.Time
  #Time_DM_min = data_idxmin.Time
  
  #Sigma_min = gb.Sigma.min()
  

  #RFIexcision.IB_Pulse_Thresh(pulses,gb,events,Sigma_min)
  #RFIexcision.Pulse_Thresh(pulses,gb,events,Sigma_min)
  
  #Clean the pulses table
  #pulses = pulses[pulses.Pulse <= RFI_percent]

  return pulses
  

def Candidates(pulses):
  pulses.groupby(['SAP','BEAM'])['Pulse'].count().groupby(level=0).mean()
  
  
  pulses.sort('Sigma',inplace=True)
  pulses.DM = pulses.DM.astype(np.float64).round(2)
  
  
  
  #Excellent
  top_count = pulses.groupby('DM')['Sigma'].count()
  top_sum = pulses.groupby('DM')['Sigma'].sum()

  cand = pd.DataFrame(columns=['DM','Sigma','N_pulses'])
  i = 0

  while not top_sum[top_sum!=0].empty:
    DM = top_sum.argmax()
    Sigma = top_sum.loc[DM-0.5:DM+0.5].sum()
    N_pulses = top_count.loc[DM-0.5:DM+0.5].count()
    cand.loc[i] = ([DM,Sigma,N_pulses])
    pulses.Candidate[(pulses.DM>=DM-0.5)&(pulses.DM<=DM-0.5)] = i
    top_count.loc[DM-0.5:DM+0.5] = 0
    top_sum.loc[DM-0.5:DM+0.5] = 0
    i += 1
    
  return 


  #Poor
  






#from scipy.optimize import curve_fit

## Define some test data which is close to Gaussian
#data = pulses.Sigma

#hist, bin_edges = numpy.histogram(data, density=True, bins=100)
#bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

## Define model function to be used to fit to the data above:
#def gauss(x, *p):
    #A, sigma = p
    #return A*numpy.exp(-x**2/(2.*sigma**2))

## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
#p0 = [7., 3.]

#coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0, maxfev=100000)

#pr = norm.ppf(1-1./(100*data.size),loc=0,scale=coeff[1])

#norm.cdf(data,loc=0,scale=coeff[2])   #IMPROVE



#from scipy.optimize import curve_fit

## Define some test data which is close to Gaussian
#data = pulses.DM

#hist, bin_edges = numpy.histogram(data, density=True, bins=1000)
#bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

## Define model function to be used to fit to the data above:
#def gauss(x, *p):
    #A, sigma = p
    #return A*numpy.exp(-x**2/(2.*sigma**2))

## p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
#p0 = [0.3, 10.]

#coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0, maxfev=100000)

#pr = norm.ppf(1-1./(100*data.size),loc=0,scale=coeff[1])








## Get the fitted curve
#hist_fit = gauss(bin_centres, *coeff)

#plt.plot(bin_centres, hist, label='Test data')
#plt.plot(bin_centres, hist_fit, label='Fitted data')

## Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
#print 'Fitted mean = ', coeff[1]
#print 'Fitted standard deviation = ', coeff[2]

#plt.show()

