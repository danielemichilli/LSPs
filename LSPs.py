#!/usr/bin/env python

'''

LOTAAS Single Pulse searcher

Written by Daniele Michilli

'''

import os
import time
import argparse

import SPclean
import LSPplot


#aggiungere verbosity e output su file e diagnostic plots

def parser():
  #---------------------
  # Command-line options
  #---------------------

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program read the database in the current directory created through DBmaker.py")
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
  #--------------------------------------------------------------------
  # If the DataBase already exists open it, otherwise creates a new one
  #--------------------------------------------------------------------
  
  args = parser()
  idL = os.path.basename(os.getcwd())
  
  #If the DataBase exists and the -w option is abstent plots it
  if (~args.w) & (os.path.isfile('SinlgePulses.hdf5')):
    print "DataBase will be plotted.\nUse -w to force the creation of a new DataBase.\n"
    LSPplot.plot(idL,args)
    
  #Otherwise create a new DataBase
  else:
    if (args.w): print "\nThe -w option has been chosen,other options will be ignored.\nA new DataBase will be created. ATTENTION: The old DataBase will be deleted!"
    else: print "\nDataBase doesn't exists"
    if 1: #raw_input("Would you create a new DataBase? It may requires many minutes. [y] n\n") == 'y':
      time0 = time.clock()  
      SPclean.obs_events(idL)  


if __name__ == '__main__':
  main()
