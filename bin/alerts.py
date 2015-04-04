import argparse
import os
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('folder', nargs=1)
args = parser.parse_args()
folder = args.folder[0]
idL = os.path.basename(os.path.dirname(folder))

if os.path.isfile('{}/ALERTS'.format(folder)):
  print 'Observation {} already analyzed!\nExiting...'.format(idL)
  exit()

pulses = pd.read_hdf('{}/SinglePulses.hdf5'.format(folder),idL+'_pulses',where=['Pulse<=1'])
pulses = pulses.loc[:,['SAP','BEAM','DM','Sigma','Pulse']]
pulses = pulses.loc[pulses.BEAM>12]
pulses = pulses.loc[pulses.DM>5]
pulses.DM = pulses.DM.astype(np.float64)

file = 0

pulses.DM = pulses.DM.round(decimals=1)
noise_level = pulses.groupby(['SAP','BEAM','DM'],sort=False).Sigma.sum().astype(np.float64)
noise_level = noise_level.median(level=['SAP','DM'])

beams = pulses.groupby(['SAP','BEAM'],sort=False)
for ind,beam in beams:
  val = beam.groupby('DM',sort=False).Sigma.sum()
  lim = noise_level.loc[ind[0]]
  lim = lim[lim.index.isin(val.index)]
  ratio = val/lim
  ratio = ratio[ratio>7]
  if not ratio.empty:
    if not file:
      file = open('{}/ALERTS'.format(folder),'w+')
      file.write('Obs. {}\n'.format(idL))
    file.write('\nPulses in SAP{}_BEAM{}\n'.format(ind[0],ind[1]))
    file.write('SAP\tBEAM\tDM\tRatio\tCounts\tSNR_max\n')
    for i in ratio.index:
      counts = beam.groupby('DM',sort=False).Sigma.count().loc[i]
      val_max = beam.groupby('DM',sort=False).Sigma.max().loc[i]
      file.write('{0:.2f}\t{1:.1f}\t{2}\t{3:.1f}\n'.format(i,ratio[i],counts,val_max))
      
  ratio = np.abs((val+5)/val.shift(10)).dropna()
  ratio = ratio[ratio>2]   
  if not ratio.empty:
    if not file:
      file = open('{}/ALERTS'.format(folder),'w+')
      file.write('Obs. {}\n\n'.format(idL))
    file.write('\nPulses in SAP{}_BEAM{}\n'.format(ind[0],ind[1]))
    file.write('SAP\tBEAM\tDM\tRatio\tCounts\tSNR_max\n')
    for i in ratio.index:
      counts = beam.groupby('DM',sort=False).Sigma.count().loc[i]
      val_max = beam.groupby('DM',sort=False).Sigma.max().loc[i]
      file.write('{0}\t{1}\t{2:.2f}\t{3:.1f}\t{4}\t{5:.1f}\n'.format(ind[0],ind[1],i,ratio[i],counts,val_max))

  

if file:
  file.close()

