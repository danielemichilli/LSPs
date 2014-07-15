import numpy as np

def SB(data):
  #Remove RFI with sinle-beam techniques
  
  print 'SB'


  DMmin = 2.
  SigmaMin = 6.
  DownfactMax = 90


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
    cond = (msk.Sigma < np.dot(msk.Sigma_inc,2)) & (np.absolute(np.subtract(msk.Time,msk.Time_inc)) < np.add(msk.Down_Time,msk.Down_Time))
    data.drop(msk.index[cond],inplace=True)

  return data


def MB(data):  #dividere tabella data in ogni sap e beam e confrontare uno alla volta
  #remove pulses in more than 1 sap
  
  print 'MB'
  data_tmp = data.copy()   #occupa molta memoria
  data_tmp['ind'] = data_tmp.index
  
  #confronta i beam anche con se stessi, aggiustare
  
  #for ind1, sap_group in data.groupby('SAP'):
  #  for ind2, event in sap_group.iterrows():
      
 
  for beam in range(13,73):
    msk = data_tmp[(data_tmp.SAP==1)|(data_tmp.SAP==2)].merge(data_tmp[(data_tmp.SAP==0) & (data_tmp.BEAM==beam)],on='DM',suffixes=['_l','_r'],copy=False)
    msk = msk[(msk.SAP_l != msk.SAP_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 2.) & (abs(msk.Time_l - msk.Time_r) < (msk.Down_Time_l + msk.Down_Time_r ))]  #provare con operazioni di numpy per aumentare efficienza
    data.drop(msk.ind_l,inplace=True)
    data.drop(msk.ind_r,inplace=True)
    
  for beam in range(13,73):
    msk = data_tmp[data_tmp.SAP==2].merge(data_tmp[(data_tmp.SAP==1) & (data_tmp.BEAM==beam)],on='DM',suffixes=['_l','_r'],copy=False)
    msk = msk[(msk.SAP_l != msk.SAP_r) & (abs(msk.Sigma_l - msk.Sigma_r) < 2.) & (abs(msk.Time_l - msk.Time_r) < (msk.Down_Time_l + msk.Down_Time_r))]  #provare con operazioni di numpy per aumentare efficienza
    data.drop(msk.ind_l,inplace=True)
    data.drop(msk.ind_r,inplace=True)
        
  return data  


def Pulses(data,grouped):
  
  data.drop(data[data.Pulse.isin(grouped[grouped.dDM>3.].index)].index,inplace=True)
  data.drop(data[data.Pulse.isin(grouped[grouped.dTime>3.*grouped.Down_Time].index)].index,inplace=True)
  
  group_temp = grouped
  
  for ind1, sap_group in group_temp[group_temp.SAP.isin([0,1])].groupby('SAP'):
    print 'SAP: ',ind1
  
    for ind2, beam_group in sap_group.groupby('BEAM'):
      
      for ind3, beam_new in group_temp[group_temp.SAP==ind1+1].groupby('BEAM'):
        
        beam_group.apply(compare,args=(beam_new,grouped),axis=1) #,raw=True)
        
      for ind3, beam_new in group_temp[group_temp.SAP==ind1+2].groupby('BEAM'):
        
        beam_group.apply(compare,args=(beam_new,grouped),axis=1) #,raw=True)        

  return data


def compare(row,beam,grouped):
  
  msk = beam[ (abs(beam.DM-row.DM)<beam.dDM+row.dDM) & (abs(beam.Time-row.Time)<beam.dTime+row.dTime) & (abs(beam.Sigma-row.Sigma)<2.) ]
    
  if len(msk)>0: 
    grouped.drop(grouped[(grouped.SAP==row.SAP) & (grouped.BEAM==row.BEAM) & (grouped.index==row.index)],inplace=True)
    grouped.drop(grouped[(grouped.SAP==beam.SAP.iloc[0]) & (grouped.BEAM==beam.BEAM.iloc[0]) & (grouped.index.isin(msk.index))],inplace=True)
  
  #  for ind2, event in sap_group.iterrows():
  
  



#AGGIUNGERE analisi dell'andamento del SNR (costante vs piccato)

#

#  return data

#usare groupby:
#In [5]: def get_letter_type(letter):
#   ...:     if letter.lower() in 'aeiou':
#   ...:         return 'vowel'
#   ...:     else:
#   ...:         return 'consonant'

#In [12]: df.groupby('A', group_keys=False).apply(lambda x: x.ix[x.B.idxmax()])
