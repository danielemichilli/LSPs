import logging
try:
  from ftplib import FTP
  import json
  import gspread
  from oauth2client.client import SignedJwtAssertionCredentials
except ImportError:
  logging.warning("Results cannot be uploaded on the website - Module missing")
import os
import time

from Paths import *


def upload(cands,folder,idL):
  upload_sheet(cands,idL)
  
  #Connect the remote server
  ftp = FTP('lsps.co.nf')
  ftp.login(user='1948849_lsps',passwd='J0140+56')

  #Create the obs folder  
  ftp.mkd(idL)

  #Create the html list file of all the obs folders
  create_list(ftp,folder)

  #Create the html obs file
  create_obs(folder)
  
  #upload the candidate plots
  ftp.cwd(idL)
  plot_list = os.listdir(folder+'/candidates/')
  for plot in plot_list:
    ftp.storbinary('STOR '+plot, open(folder+'/candidates/'+plot,'r'))
  
  ftp.quit()
  os.remove(folder+'/obs_list.html')
  os.remove(folder+'candidates/observation.html')
  
  return


def create_list(ftp,folder):
  obs_list = []
  ftp.retrlines('MLSD',obs_list.append)
  obs_list = [line[1:] for line in obs_list if (line[1]=='L')&(line[2:].isdigit())]
  obs_list.sort()  
  
  f = open(folder+'/obs_list.html','w')
  f.write('<!DOCTYPE html>\n<html>\n<head>\n<base target="plots">\n</head>\n<body>\n <ul style="list-style-type:none">\n\n')
  for line in obs_list: 
    s = '  <li><a href="{idL}/observation.html" target="plots">{idL}</a></li>'.format(idL=line)
    f.write(s)
  f.write('\n\n </ul>\n</body>\n</html>\n')
  f.close()
  
  ftp.storbinary('STOR obs_list.html', open(folder+'/obs_list.html','r'))
  
  return


def create_obs(folder):
  f = open(folder+'candidates/observation.html','w')
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
    
  plot_list = os.listdir(folder+'/candidates/')
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
  
  for idx,cand in cands.iterrows():
    if cand.N_pulses == 1: kind = 'SP'
    else: kind = 'RC'
    row = [cands.id, date, idL, cands.SAP, cands.BEAM, kind, cands.N_pulses, cands.DM, cands.Rank]
    wks.append_row(row)
  
  return
