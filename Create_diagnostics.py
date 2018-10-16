import argparse
import os
import shutil

import pandas as pd

import src.Paths as PATH
from LSpS import set_paths


def parser():
  '''
  Command-line options
  '''

  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program creates a DataBase of the events in the observation.")
  parser.add_argument('id_obs', help='Observation ID')
  parser.add_argument('-plot', help="Run over confirmation observations.",action='store_true')
  parser.add_argument('-store_online', help="Run over confirmation observations.",action='store_true')
  parser.add_argument('-conf', help="Run over confirmation observations.",action='store_true')
  args = parser.parse_args()
  
  return args
  

def load_DB():
  meta_data = pd.read_hdf(PATH.DB, 'meta_data')
  pulses = pd.read_hdf(PATH.DB, 'pulses')
  cands = pd.read_hdf(PATH.DB, 'candidates')
  cands = cands[cands.main_cand == 0]
  cands.sort_values('Sigma', inplace=True, ascending=False)
  cands = cands.groupby('BEAM').head(10)
  cands = cands.head(50)
  cands = cands[ ((cands.N_pulses == 1) & (cands.Sigma>10.)) | ((cands.N_pulses > 1) & (cands.Sigma>16.)) ]
  cands.sort_values('Sigma', inplace=True, ascending=False)
  return 


def main(PATH):
  args = parser()
  PATH = set_paths(args, PATH)
  PATH.DB = os.path.join(PATH.OBS_FOLDER,'sp/SinglePulses.hdf5')
  if os.path.isdir(os.path.join(PATH.OBS_FOLDER, 'sp/candidates')): shutil.rmtree(os.path.join(PATH.OBS_FOLDER, 'sp/candidates'))
  
  
  if args.conf: inc = 0
  else: inc = 12

  meta_data, pulses, cands = load_DB()
  if args.plot: LSPplot.output(args.id_obs, pulses, meta_data, cands, PATH.DB, inc=inc)
  if args.store_online: Internet.upload(cands, args.id_obs, os.path.join(PATH.WRK_FOLDER,'sp/candidates/.'), meta_data, pulses)
  return


if __name__ == '__main__':
  try: main(PATH)
  finally:
    shutil.copytree(os.path.join(PATH.WRK_FOLDER, 'sp/candidates'), os.path.join(PATH.OBS_FOLDER, 'sp/candidates'))
    shutil.rmtree(PATH.WRK_FOLDER)
    shutil.rmtree(PATH.TMP_FOLDER)
