import logging
try:
  from ftplib import FTP
  from ftplib import error_perm
  import json
  import gspread
  from oauth2client.client import SignedJwtAssertionCredentials
except ImportError:
  logging.warning("Results cannot be uploaded on the website - Module missing")
import os
import time

from Paths import *


def upload(cands,folder,idL):
  #Retrive website information to connect
  f = open(SITE_ATTR)
  site_host = f.readline().split()[0]
  username = f.readline().split()[0]
  psw = f.readline().split()[0]
  f.close()
  
  #Connect the remote server
  ftp = FTP(site_host)
  ftp.login(user=username,passwd=psw)
  ftp.cwd('public_html')
  
  #Create the obs folder  
  try: ftp.mkd(idL)
  except error_perm:
    ftp.cwd(idL)
    for file in ftp.nlst():
      try: ftp.delete(file)
      except error_perm: pass
    ftp.cwd('..')
    ftp.rmd(idL)
    ftp.mkd(idL)
  
  #Create the html list file of all the obs folders
  create_list(ftp,folder,idL)

  #Create the html obs file
  create_obs(folder,idL)
  
  #upload the candidate plots
  ftp.cwd(idL)
  plot_list = os.listdir(folder+idL+'/sp/candidates/')
  for plot in plot_list:
    ftp.storbinary('STOR '+plot, open(folder+idL+'/sp/candidates/'+plot,'r'))
  
  ftp.quit()
  os.remove(folder+idL+'/sp/obs_list.html')
  os.remove(folder+idL+'/sp/candidates/observation.html')
  
  return


def create_list(ftp,folder,idL):
  #obs_list = []
  #ftp.retrlines('MLSD',obs_list.append)
  obs_list = ftp.nlst()
  obs_list = [line for line in obs_list if (line[0]=='L')&(line[2:].isdigit())]
  obs_list.sort()  
  
  f = open(folder+idL+'/sp/obs_list.html','w')
  f.write('<!DOCTYPE html>\n<html>\n<head>\n<base target="plots">\n</head>\n<body>\n <ul style="list-style-type:none">\n\n')
  for line in obs_list: 
    s = '  <li><a href="{idL}/observation.html" target="plots">{idL}</a></li>'.format(idL=line)
    f.write(s)
  f.write('\n\n </ul>\n</body>\n</html>\n')
  f.close()
  
  ftp.storbinary('STOR obs_list.html', open(folder+idL+'/sp/obs_list.html','r'))
  
  return


def create_obs(folder,idL):
  f = open(folder+idL+'/sp/candidates/observation.html','w')
  f.write("""<!DOCTYPE html>
<html>
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  <script type="text/javascript">
  //<![CDATA[
  try{if (!window.CloudFlare) {var CloudFlare=[{verbose:0,p:0,byc:0,owlid:"cf",bag2:1,mirage2:0,oracle:0,paths:{cloudflare:"/cdn-cgi/nexp/dok3v=1613a3a185/"},atok:"5ccac93e3287bbdceb1966dca9ac87fc",petok:"97f84cf4d5bb7697926ef4eb5aac680b571c8f73-1441792405-1800",betok:"0e2c131e1500db8ab0fdb5c4f6d43a9b63524ce1-1441792405-120",zone:"hawksey.info",rocket:"a",apps:{}}];document.write('<script type="text/javascript" src="//ajax.cloudflare.com/cdn-cgi/nexp/dok3v=e9627cd26a/cloudflare.min.js"><'+'\/script>');}}catch(e){};
  //]]>
  </script>
  <script data-cfasync="false" src="//ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
  <script data-cfasync="false" type="text/javascript" src="../send_to_sheet.js"></script>
</head>
<body>

""")
    
  plot_list = os.listdir(folder+idL+'/sp/candidates/')
  for plot in plot_list:
    if os.path.splitext(plot)[1] == ".png":
      f.write('<form id="foo">\n')
      f.write(' <input type="text" name="id" value="{}" readonly>\n'.format(os.path.splitext(plot)[0]))
      f.write(' <select name="type">\n  <option value="rfi">RFI</option>\n  <option value="poor">Poor candidate</option>\n  <option value="excellent">Excellent candidate</option>\n  <option value="known">Known source</option>\n </select>\n')
      f.write(' <br>\n Notes: <input type="text" name="notes">\n <br><br>\n <input type="submit" value="Send"/>\n </form>\n')
      f.write('<img src="{}" width="100%">\n'.format(plot))
      f.write('<hr>\n\n')
  
  f.write('\n\n<script data-cfasync="false" type="text/javascript" src="../send_to_sheet.js"></script>\n</body>\n</html>')
  f.close()
  
  return


def upload_sheet(cands,idL):
  try: json_key = json.load(open(SITE_CERT))
  except IOError:
    logging.warning("Spreadsheet cannot be uploaded - Google certificate missing")
    return
  scope = ['https://spreadsheets.google.com/feeds']
  credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
  gc = gspread.authorize(credentials)
  wks = gc.open("LSPs Candidates").sheet1
  date = time.strftime("%d/%m/%Y")
  
  cand_list = wks.col_values(1)
  
  for idx,cand in cands.iterrows():
    if cand.N_pulses == 1: kind = 'SP'
    else: kind = 'RC'
    row = [cand.id, date, idL, cand.SAP, cand.BEAM, kind, cand.N_pulses, cand.DM, cand.Rank]
    if cand.id in cand_list:    #Magari meglio rimuovere tutta l'osservazione nel sheet
      row_num = wks.find(cand.id).row
      col_num = len(wks.row_values(row_num))
      #col_num = chr(len(row)-1 + ord('A'))
      col_num = chr(col_num-1 + ord('A'))
      cells = wks.range('A{row}:{col}{row}'.format(row=row_num,col=col_num))
      for i,cell in enumerate(cells): 
        if i < len(row): cell.value = row[i]
        else:  cell.value = ''
      wks.update_cells(cells)
      logging.warning("Candidate {} already existed in the datasheet! It has been substituted.".format(cand.id))
    else:
      wks.append_row(row)
  
  return
