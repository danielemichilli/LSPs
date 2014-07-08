#!/usr/bin/env python

'''
Read data from singlepulse beam files, transfer it on the memory, elaborate it and store at the end.

Written by Daniele Michilli
'''


import os
import time

import SPclean
from ssps import ssps_grab_lotaas



def main():
  
  idL = os.path.basename(os.getcwd())
  
  #Write the empty DB if it doesn't exist
  #if 'SinlgePulses.hdf5': 
  #  print "The DataBase already exists."
  #  exit()
    
  data = SPclean.obs_events(idL)
  
  data.sort(['SAP','BEAM','DM','Time','Sigma'],inplace=True)
  
  
  for sap in range(0,3):  #0,3
    for beam in range(13,74):  #13,74
      SB = data[(data.SAP==sap)&(data.BEAM==beam)]
      if not SB.empty:
        ssps_grab_lotaas.ssps_grab(idL,data[(data.SAP==sap)&(data.BEAM==beam)])
  
  
  



if __name__ == '__main__':
  main()
