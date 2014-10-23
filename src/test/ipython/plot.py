#######################################################
#
# Written by Daniele Michilli
#
########################################################

import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import numpy as np

#--------------------
#Command-line options
#--------------------

#def plot():

idL = os.path.basename(os.getcwd())



try:
  store = pd.HDFStore('SinlgePulses.hdf5','r')
except IOError:
  print "DataBase doesn't exist in the current folder."
  exit()
grouped = store[idL+'_pulses']
data = store[idL]

cond1 = (grouped.DM/(grouped.DM_min+grouped.dDM) > 0.095) & (grouped.DM/(grouped.DM_min+grouped.dDM) < 1.005)
cond2 = grouped.Sigma/grouped.Sigma_min> 1.1
#grouped = grouped[cond1 & cond2]
#pulse = grouped[(grouped.DM>43) & (grouped.DM<44)]
#rfi = grouped[(grouped.DM<10)]
#grouped = grouped[(grouped['DM']>dmlo) & (grouped['DM']<dmhi) & (grouped['Time']>tmlo) & (grouped['Time']<tmhi) & (grouped['SAP'].isin(sap)) & (grouped['BEAM'].isin(beam))]
#pulse['ratio'] = pulse.Sigma/pulse.Sigma_min
#rfi['ratio'] = rfi.Sigma/rfi.Sigma_min
#grouped.sort('ratio')
#print 'Ratio: ',float(len(rfi))/float(len(pulse))
#print 'Tot: ',len(pulse)


grouped1 = grouped[grouped.DM<10]
grouped2 = grouped[(grouped.DM>43.46)&(grouped.DM<43.54)]
grouped3 = grouped[grouped.DM>50]


#grouped1['dm_shift'] = grouped1.DM/(grouped1.DM_min+grouped1.dDM)
#grouped2['dm_shift'] = grouped2.DM/(grouped2.DM_min+grouped2.dDM)
#grouped3['dm_shift'] = grouped3.DM/(grouped3.DM_min+grouped3.dDM)

#grouped1['dm_shift'] = grouped1.Sigma/grouped1.Sigma_min
#grouped2['dm_shift'] = grouped2.Sigma/grouped2.Sigma_min
#grouped3['dm_shift'] = grouped3.Sigma/grouped3.Sigma_min

grouped1['dm_shift'] = grouped1.dTime/grouped1.dDM
grouped2['dm_shift'] = grouped2.dTime/grouped2.dDM
grouped3['dm_shift'] = grouped3.dTime/grouped3.dDM

#grouped1['dm_shift'] = grouped1.Sigma_DM_max/grouped1.Sigma_DM_min
#grouped2['dm_shift'] = grouped2.Sigma_DM_max/grouped2.Sigma_DM_min
#grouped3['dm_shift'] = grouped3.Sigma_DM_max/grouped3.Sigma_DM_min

#grouped1['dm_shift'] = grouped1.Sigma/grouped1.Sigma_DM_min
#grouped2['dm_shift'] = grouped2.Sigma/grouped2.Sigma_DM_min
#grouped3['dm_shift'] = grouped3.Sigma/grouped3.Sigma_DM_min


lim_min = 0
lim_max = 7

first1 = len(grouped1)
first2 = len(grouped2)
first3 = len(grouped3)

grouped1 = grouped1[(grouped1.dm_shift<lim_max)&(grouped1.dm_shift>lim_min)]
grouped2 = grouped2[(grouped2.dm_shift<lim_max)&(grouped2.dm_shift>lim_min)]
grouped3 = grouped3[(grouped3.dm_shift<lim_max)&(grouped3.dm_shift>lim_min)]

print grouped2[(grouped2.dm_shift>=lim_max)|(grouped2.dm_shift<=lim_min)]

print ' Low-DM RFI survived: ', len(grouped1)/float(first1)*100.
print 'High-DM RFI survived: ', len(grouped3)/float(first3)*100.
print 'Real pulses survived: ', len(grouped2)/float(first2)*100.

#grouped1['dm_shift'].hist(bins=200,alpha=.5)
grouped2['dm_shift'].hist(bins=200,alpha=.5)
#grouped3['dm_shift'].hist(bins=200,alpha=.5)

#plt.plot(grouped.Time,grouped.DM,'k.')  


plt.xlabel('DM')
plt.ylabel('Ratio')

plt.show()


store.close()