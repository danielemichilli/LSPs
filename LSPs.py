#!/usr/bin/env python

'''
Read data from singlepulse beam files, transfer it on the memory, elaborate it and store at the end.

Written by Daniele Michilli
'''


import os
import time

import SPclean
from ssps import ssps_grab_lotaas
import Group


def main():
  
  idL = os.path.basename(os.getcwd())
  
  #-------------------------------------------------------------------
  #If the DataBase already exists open it, otherwise creates a new one
  #-------------------------------------------------------------------
  
  #if os.path.isfile('SinlgePulses.hdf5'):
  #  print "DataBase already exists in the current folder, it will be plotted.\nUse -w to force the creation of a new DataBase.\n"
  #  store = pd.HDFStore('SinlgePulses.hdf5','r')
  #  data = store[idL]
  #  store.close()
  #else:
  #  
    
  time0 =time.clock()  
  data = SPclean.obs_events(idL)
  print 'Time: ',time.clock()-time0,' s'
  
#  for sap in range(0,3):  #0,3
#    for beam in range(13,74):  #13,74
#      SB = data[(data.SAP==sap)&(data.BEAM==beam)]
#      if not SB.empty:
#        ssps_grab_lotaas.ssps_grab(idL,data[(data.SAP==sap)&(data.BEAM==beam)])
  
  
  



if __name__ == '__main__':
  main()
