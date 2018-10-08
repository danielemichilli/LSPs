import os
import time
import argparse
import shutil
import subprocess
import warnings

import numpy as np
import pandas as pd

from src import SPclean
from src.Paths import *


def parser():
  '''
  Command-line options
  '''

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program creates a DataBase of the events in the observation.")
  parser.add_argument('id_obs', help='Observation ID')
  parser.add_argument('-conf', help="Run over confirmation observations.",action='store_true')
  parser.add_argument('-load_DB', help="Load events from a DataBase instead of text files.",action='store_true')
  parser.add_argument('--debug', help="Debug mode: shows all the warnings.",action='store_true')
  args = parser.parse_args()
  
  return args


def main():
  #warnings.filterwarnings('error', category=FutureWarning)

  args = parser()
  if args.conf:
    global OBS_FOLDER
    global RAW_FOLDER
    OBS_FOLDER = os.path.join(OBS_FOLDER, 'confirmations/')
    RAW_FOLDER = os.path.join(RAW_FOLDER, 'confirmations/')
  args.folder = OBS_FOLDER

  try: shutil.rmtree(WRK_FOLDER.format(args.id_obs))
  except OSError: pass
  os.makedirs(WRK_FOLDER.format(args.id_obs))
  os.makedirs('{}/sp'.format(WRK_FOLDER.format(args.id_obs)))
  os.makedirs('{}/sp/candidates'.format(WRK_FOLDER.format(args.id_obs)))
  try: shutil.rmtree(TMP_FOLDER.format(args.id_obs))
  except OSError: pass
  os.makedirs(TMP_FOLDER.format(args.id_obs))
  try: shutil.rmtree('{}{}/sp'.format(OBS_FOLDER,args.id_obs))
  except OSError: pass

  scriptFolder = os.path.dirname(os.path.realpath(__file__))
  git_folder = os.path.join(scriptFolder, '.git')
  vers = subprocess.check_output(['git','--git-dir',git_folder,'describe','--tags','--abbrev=0','--always']).strip()
  SPclean.log("L-SpS version used: {}".format(vers), args.id_obs)
  args.vers = vers

  time0 = time.time()
  SPclean.log("The DataBase is being created", args.id_obs)

  try:
    SPclean.main(args)
    SPclean.log("The DataBase has been created", args.id_obs)
    SPclean.log("Time spent: {:.2f} s".format(time.time() - time0), args.id_obs)
  
  except:
    SPclean.log_err("Fatal: an error arised in processing the observation", args.id_obs)

  finally:
    shutil.copytree('{}/sp'.format(WRK_FOLDER.format(args.id_obs)),'{}/{}/sp'.format(OBS_FOLDER,args.id_obs))
    shutil.rmtree(WRK_FOLDER.format(args.id_obs))
    shutil.rmtree(TMP_FOLDER.format(args.id_obs))

  return


if __name__ == '__main__':
  main()
