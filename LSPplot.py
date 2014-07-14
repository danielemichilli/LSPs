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

idL = os.path.basename(os.getcwd())

#Define the command-line parsers
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="The program read the database in the current directory created through DBmaker.py")
#parser.add_argument("-p", type=int, nargs=1, default=0, help="Plot level: 0-no plot, 1-plot on screen, 2-plot on file, 3-plot on screen AND file.")
#parser.add_argument("-v", help="Print the output on screen.")
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

#Assign names to the parsers
dmlo = args.dmlo[0]
dmhi = args.dmhi[0]
tmlo = args.tmlo[0]
tmhi = args.tmhi[0]
sap = args.sap
beam = args.beam


#-----------------------------------
#Selection of values in the database
#-----------------------------------

try:
  store = pd.HDFStore('SinlgePulses.hdf5','r')
except IOError:
  print "DataBase doesn't exist in the current folder."
  exit()

data = store[idL]

print len(data)

#provare se e piu veloce (mettere default=[-1] nei parser):
#if dmlo>0:
#  data=data[(data['DM']>dmlo)]
#if dmhi>0:
#  data=data[(data['DM']<dmhi)]
#...

data=data[(data['DM']>dmlo) & (data['DM']<dmhi) & (data['Time']>tmlo) & (data['Time']<tmhi) & (data['SAP'].isin(sap)) & (data['BEAM'].isin(beam))]

if args.s: 
  sig=(data.Sigma/5.)**3
else:
  sig=8

if args.l: 
  linewidths=data.Downfact*0.1
else:
  linewidths=1
  
if args.c: 
  #col = data.SAP *10 + (data.BEAM-13) *10./62.
  
  col = data.Pulse.astype(float)
  
else:
  col='b'  
  
print data
  
plt.scatter(data.Time, data.DM, facecolors='none', s=sig, c=col, cmap=mpl.cm.rainbow)

#if args.c:   #Testare che faccia tutto bene, sembra troppo robusto
#  ticks = np.linspace(col.min(),col.max(),num=10)
#  bar = plt.colorbar(ticks=ticks)
#  bar.set_ticklabels(['{0:.0f}, {1:.0f}'.format(int(t)/10,t%10/10.*62.+13) for t in ticks])
#  bar.ax.set_xlabel('sap, beam',ha='left',labelpad=-380)
#  bar.update_ticks
#  bar.ax.xaxis.set_ticks_position('top')
  
  
#confrontare plot e scatter: velocita e bellezza
#plt.plot(data['Time'], data['DM'], 'ko', mfc='none', ms=2)

plt.xlabel('Time (s)')
plt.ylabel('DM (pc/cm3)')
plt.axis([tmlo,tmhi,dmlo,dmhi])

plt.show()


store.close()