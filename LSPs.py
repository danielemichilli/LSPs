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
  
  
  #grab(idL
  
  
  



if __name__ == '__main__':
    main()
