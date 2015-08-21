import os

import LSPplot


def output(folder,idL,pulses,events,meta_data):
  pulses.sort('Sigma',ascending=False,inplace=True)
  
  store = '{}{}/sp/beam_plots'.format(folder,idL)
  
  gb_puls = pulses.groupby(['SAP','BEAM'],sort=False)
  
  for n in gb_puls.indices.iterkeys():
    name = 'SAP{}_BEAM{}'.format(n[0],n[1])
    os.makedirs('{}/{}'.format(store,name))
    
  #pool = mp.Pool()
  #pool.map(output_beams, [(pulses,events,meta_data,store,idL,n) for n in gb_puls.indices.iterkeys()])
  #pool.close()
  #pool.join()  
  
  for n in gb_puls.indices.iterkeys():
    output_beams((pulses,events,meta_data,store,idL,n))
  
  pulses = pulses[pulses.Pulse==0]
  output_pointing(pulses,folder,idL)
  
  return



def output_beams((pulses,events,meta_data,folder,obs,(sap,beam))):
  
  top = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==0)]
  good = pulses[(pulses.SAP==sap) & (pulses.BEAM==beam) & (pulses.Pulse==1)]
  
  if beam == 12:
    if top.shape[0] > 0: 
      LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates(0-9).png'.format(folder,sap,beam),obs)
      if top.shape[0] > 10:
        LSPplot.sp_shape(top.iloc[10:20],events,'{}/SAP{}_BEAM{}/top_candidates(10-19).png'.format(folder,sap,beam),obs)
        if top.shape[0] > 20: 
          LSPplot.sp_shape(top.iloc[20:30],events,'{}/SAP{}_BEAM{}/top_candidates(20-29).png'.format(folder,sap,beam),obs)
    LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
    
  else:
    if top.shape[0] > 0: 
      LSPplot.sp_shape(top.head(10),events,'{}/SAP{}_BEAM{}/top_candidates.png'.format(folder,sap,beam),obs)
      LSPplot.sp_plot(top,good,meta_data,sap,beam,'{}/SAP{}_BEAM{}/beam.png'.format(folder,sap,beam))
  
  return
  
  

def output_pointing(pulses,folder,idL):
  
  top_candidates = pulses[pulses.BEAM>12].groupby(['SAP','BEAM'],sort=False).head(10)
  top_candidates = top_candidates.append(pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30),ignore_index=False)
  top_candidates.sort(['SAP','BEAM','Sigma'],ascending=[True,True,False],inplace=True)
  top_candidates['code'] = top_candidates.index
  if not top_candidates.empty:
    a = top_candidates.groupby(['SAP','BEAM'],sort=False).apply(lambda x: range(len(x))).tolist()
    b = [val for sublist in a for val in sublist]
    top_candidates.index = b
  top_candidates.Duration *= 1000
  top_candidates['void'] = ''
  top_candidates.to_csv('{}{}/sp/files/top_candidates.inf'.format(folder,idL),sep='\t',float_format='%.2f',\
    columns=['SAP','BEAM','Sigma','DM','void','Time','void','Duration','void','code'],\
    header=['SAP','BEAM','Sigma','DM (pc/cm3)','Time (s)','Duration (ms)','code','','',''],index_label='rank',encoding='utf-8')
  
  puls = pulses[pulses.BEAM==12].groupby('SAP',sort=False).head(30)
  if not puls.empty: LSPplot.obs_top_candidates(puls,store='{}{}/sp/files/inc_top_candidates.png'.format(folder,idL),incoherent=True)
  
  for sap in pulses.SAP.unique():
    puls = pulses[(pulses.SAP==sap)&(pulses.BEAM>12)].groupby('BEAM',sort=False).head(10)
    if not puls.empty: LSPplot.obs_top_candidates(puls,store='{}{}/sp/files/top_candidates(SAP{}).png'.format(folder,idL,sap))

  return


 
