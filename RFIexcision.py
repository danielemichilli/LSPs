import numpy as np

DMmin = 3.
SigmaMin = 6.
DownfactMax = 90


def SB(data):
  #Remove RFI with sinle-beam techniques
  
  print 'SB'

  #Remove low DM
  data = data[data.DM>DMmin]
  
  #Remove low sigma
  data = data[data.Sigma>SigmaMin]

  #Remove high downfactors
  data = data[data.Downfact<DownfactMax]
  
  return data


def IB(data,incoher):  #confrontare anche dm vicine
  #Remove RFI with multi-beam techniques
  
  print 'IB'
  
  #incoherent beam
  if not incoher.empty:
    
    msk = data.merge(incoher,on='DM',suffixes=['','_inc'],copy=False,right_index=True)
    cond = (msk.Sigma < np.dot(msk.Sigma_inc,2)) & (np.absolute(np.subtract(msk.Time,msk.Time_inc)) < np.multiply(np.add(np.multiply(msk.Downfact,msk.Sampling),np.multiply(msk.Downfact_inc,msk.Sampling_inc)),2))
    data.drop(msk.index[cond],inplace=True)

  return data


def MB(data):  #dividere tabella data in ogni sap e beam e confrontare uno alla volta
  #remove pulses in more than 1 sap
  
  print 'MB'
  data_tmp = data.copy()   #occupa molta memoria
  data_tmp['ind'] = data_tmp.index
  
  #confronta i beam anche con se stessi, aggiustare
  
 
  for beam in range(13,73):
    msk = data_tmp[(data_tmp.SAP==1)|(data_tmp.SAP==2)].merge(data_tmp[(data_tmp.SAP==0) & (data_tmp.BEAM==beam)],on='DM',suffixes=['_l','_r'],copy=False)
    msk = msk[(msk.SAP_l != msk.SAP_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 1.) & (abs(msk.Time_l - msk.Time_r) < (2. * (msk.Downfact_l * msk.Sampling_l + msk.Downfact_r * msk.Sampling_r)))]  #provare con operazioni di numpy per aumentare efficienza
    #msk = msk[(msk.BEAM_l != msk.BEAM_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 3.) & ((abs(msk.Time_l - msk.Time_r)) < (2. * (msk.Downfact_l * msk.Sampling_l + msk.Downfact_r * msk.Sampling_r)))]  #SOLO PER I TEST!
    data.drop(msk.ind_l,inplace=True)
    data.drop(msk.ind_r,inplace=True)
    
  for beam in range(13,73):
    msk = data_tmp[data_tmp.SAP==2].merge(data_tmp[(data_tmp.SAP==1) & (data_tmp.BEAM==beam)],on='DM',suffixes=['_l','_r'],copy=False)
    msk = msk[(msk.SAP_l != msk.SAP_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 2.) & (abs(msk.Time_l - msk.Time_r) < (2. * (msk.Downfact_l * msk.Sampling_l + msk.Downfact_r * msk.Sampling_r)))]  #provare con operazioni di numpy per aumentare efficienza
    #msk = msk[(msk.BEAM_l != msk.BEAM_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 3.) & ((abs(msk.Time_l - msk.Time_r)) < (2. * (msk.Downfact_l * msk.Sampling_l + msk.Downfact_r * msk.Sampling_r)))]  #SOLO PER I TEST!
    data.drop(msk.ind_l,inplace=True)
    data.drop(msk.ind_r,inplace=True)
        
  return data  


#def Isolated(data):   #forse prima sort in DM and Time, poi slice up and down 
#  #Remove isolate events
#  msk = data
#  msk.DM[msk.DM<=40.] = np.trunc(msk.DM[msk.DM<=40.],0)  #mask,truncate
#  msk.Time = np.trunc(msk.Time,3)
#  msk = msk.groupby([DM,Time],sort=false).count > 1
#  data.drop(msk.index,inplace=True)
#  
#  data.sort(['DM','Time','Sigma'],inplace=True)
#  data['DM'<40].drop((data.diff()>0.3) & (data['Time'].diff()>3.*DownfactLow) )
#
#
#def myround(x, base=5):
#    return int(base * round(float(x)/base))
  
  


#

#  return data

#usare groupby:
#In [5]: def get_letter_type(letter):
#   ...:     if letter.lower() in 'aeiou':
#   ...:         return 'vowel'
#   ...:     else:
#   ...:         return 'consonant'

#In [12]: df.groupby('A', group_keys=False).apply(lambda x: x.ix[x.B.idxmax()])
