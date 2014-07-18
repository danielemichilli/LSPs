#!/usr/bin/env python

'''
Read data from singlepulse beam files, transfer it on the memory, elaborate it and store at the end.

Written by Daniele Michilli
'''


import os
import time
import argparse

import SPclean
import LSPplot


def parser():
  #--------------------
  #Command-line options
  #--------------------

  #Define the command-line parsers
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program read the database in the current directory created through DBmaker.py")
  #parser.add_argument("-p", type=int, nargs=1, default=0, help="Plot level: 0-no plot, 1-plot on screen, 2-plot on file, 3-plot on screen AND file.")
  #parser.add_argument("-v", help="Print the output on screen.")
  parser.add_argument("-w", help="Force the creation of a DataBase also if another is already present.", action='store_true')
  parser.add_argument("-s", help="Plot marks with size depending on sigma (slower).", action='store_true')
  parser.add_argument("-l", help="Plot marks with line width depending on down factor (slower).", action='store_true')
  parser.add_argument("-c", help="Plot marks with colour depending on beam number (slower).", action='store_true')
  parser.add_argument("-dmlo", type=float, nargs=1, default=[0.0], help="Set the lower DM to show.")
  parser.add_argument("-dmhi", type=float, nargs=1, default=[500.0], help="Set the higher DM to show.")
  parser.add_argument("-tmlo", type=float, nargs=1, default=[0.0], help="Set the lower time to show.")
  parser.add_argument("-tmhi", type=float, nargs=1, default=[3600.0], help="Set the higher time to show.")
  parser.add_argument("-sap", type=int, nargs='+', default=range(0,3), help="Set the SAP number.")
  parser.add_argument("-beam", type=int, nargs='+', default=range(13,74), help="Set the BEAM number.")
  args = parser.parse_args()

  return args




def main():
  
  args = parser()
  idL = os.path.basename(os.getcwd())
  
  #--------------------------------------------------------------------
  # If the DataBase already exists open it, otherwise creates a new one
  #--------------------------------------------------------------------
  
  if (~args.w) & (os.path.isfile('SinlgePulses.hdf5')):
    print "DataBase will be plotted.\nUse -w to force the creation of a new DataBase.\n"
    LSPplot.plot(idL,args)
    
  
  else:
    if (args.w): print "\nThe -w option has been chosen,other options will be ignored.\nA new DataBase will be created. ATTENTION: The old DataBase will be deleted!"
    else: print "\nDataBase doesn't exists"
    if raw_input("Would you create a new DataBase? It may requires many minutes. [y] n\n") == 'y':
      time0 = time.clock()  
      data = SPclean.obs_events(idL)
      print 'Time: ',time.clock()-time0,' s'
  


if __name__ == '__main__':
  main()
