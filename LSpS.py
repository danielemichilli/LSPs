import os
import time
import argparse
import shutil
import subprocess
import warnings
import sys
import traceback
from StringIO import StringIO

import numpy as np
import pandas as pd

from src import SPclean
import src.Paths as PATH


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


def main(PATH):
  #warnings.filterwarnings('error', category=FutureWarning)

  args = parser()
  if args.conf:
    PATH.OBS_FOLDER = os.path.join(PATH.OBS_FOLDER, 'confirmations')
    PATH.RAW_FOLDER = os.path.join(PATH.RAW_FOLDER, 'confirmations')
  
  PATH.WRK_FOLDER = os.path.join(PATH.WRK_FOLDER, args.id_obs)
  PATH.TMP_FOLDER = os.path.join(PATH.TMP_FOLDER, args.id_obs)
  PATH.OBS_FOLDER = os.path.join(PATH.OBS_FOLDER, args.id_obs)

  if os.path.isdir(PATH.WRK_FOLDER): shutil.rmtree(PATH.WRK_FOLDER)
  os.makedirs(PATH.WRK_FOLDER)
  os.makedirs(os.path.join(PATH.WRK_FOLDER, 'sp'))
  os.makedirs(os.path.join(PATH.WRK_FOLDER, 'sp/candidates'))
  if os.path.isdir(PATH.TMP_FOLDER): shutil.rmtree(PATH.TMP_FOLDER)
  os.makedirs(PATH.TMP_FOLDER)
  if os.path.isdir(os.path.join(PATH.OBS_FOLDER, 'sp')): shutil.rmtree(os.path.join(PATH.OBS_FOLDER, 'sp'))

  result = StringIO()
  sys.stdout = result

  scriptFolder = os.path.dirname(os.path.realpath(__file__))
  git_folder = os.path.join(scriptFolder, '.git')
  vers = subprocess.check_output(['git','--git-dir',git_folder,'describe','--tags','--abbrev=0','--always']).strip()
  print "L-SpS version used: {}".format(vers)
  args.vers = vers

  time0 = time.time()
  print "The DataBase is being created."

  try:
    SPclean.main(args)
    print "The DataBase has been created."
    print "Time spent: {:.2f} s.".format(time.time() - time0)
    
  except Exception:
    with open(os.path.join(PATH.WRK_FOLDER, 'sp/ERROR_log.txt'), 'w') as stderr:
      stderr.write(traceback.format_exc())

  finally:
    shutil.copytree(os.path.join(PATH.WRK_FOLDER, 'sp'), os.path.join(PATH.OBS_FOLDER, 'sp'))
    with open(os.path.join(PATH.OBS_FOLDER, 'sp/log.txt'), 'w') as log:
      log.write(result.getvalue())
    shutil.rmtree(PATH.WRK_FOLDER)
    shutil.rmtree(PATH.TMP_FOLDER)
    stdout.close()
    stderr.close()

  return


if __name__ == '__main__':
  main(PATH)
