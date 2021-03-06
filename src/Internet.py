from requests import ConnectionError
import logging
import subprocess
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import os
import time

import Parameters
import Paths as PATH


def upload(cands,idL,folder,meta_data,pulses):
  upload_sheet(cands,idL,meta_data,pulses)
  upload_plots(idL,folder)
  return


def upload_plots(idL,folder):
  FNULL = open(os.devnull, 'w')
  error = subprocess.call(['scp','-prq',folder,'ag004:/var/www/lofarpwg/lotaas-sp/observations/{}/'.format(idL)], stdout=FNULL, stderr=FNULL)
  if error: 
    subprocess.call(['ssh','ag004','rm -rf /var/www/lofarpwg/lotaas-sp/observations/{}'.format(idL)], stdout=FNULL, stderr=FNULL)
    raise ConnectionError("ATTENTION! Impossible to upload plots to the website. Try to upload the observation later with Upload.py")
  return


def upload_sheet(cands,idL,meta_data,pulses):
  scope = ['https://spreadsheets.google.com/feeds']
  certificate = os.path.join(PATH.PL_FOLDER, "LSP/.site_cert.json")
  credentials = ServiceAccountCredentials.from_json_keyfile_name(certificate, scope)
  gc = gspread.authorize(credentials)
  sh = gc.open("LSP_candidates")
  wks = sh.worksheet("Classifier_v4")
  col_size = wks.col_count
  col_size_letter = chr(col_size-1 + ord('A'))
  row_size = wks.row_count

  #Remove old rows for the same observation
  obs_list = wks.col_values(4)
  row_list = [i+1 for i,val in enumerate(obs_list) if val == idL]
  for row_num in row_list:
    cells = wks.range('A{row}:{col}{row}'.format(row=row_num,col=col_size_letter))
    for cell in cells : cell.value = ''
    wks.update_cells(cells)  #Possible to append all together (faster)?
      
  #Append candidates
  date = time.strftime("%m/%d/%Y")
  fileName, fileExtension = os.path.splitext(Parameters.__file__)
  git_folder = '{}/.git'.format(os.path.dirname(os.path.dirname(fileName)))
  vers = subprocess.check_output(['git','--git-dir',git_folder,'describe','--tags','--abbrev=0','--always']).strip()
  for idx,cand in cands.iterrows():
    meta = meta_data[(meta_data.SAP == cand.SAP) & (meta_data.BEAM == cand.BEAM)]
    link = '=HYPERLINK(CONCATENATE("http://www.astron.nl/lofarpwg/lotaas-sp/observations/{}/";OFFSET($A$1;ROW()-1;0);".pdf");"Plot")'.format(idL)
    RA = "'" + meta.RA[0].split('.')[0]
    DEC = "'" + meta.DEC[0].split('.')[0]
    cand_percent = int(float(cand.N_pulses) / pulses[(pulses.Pulse==0) & (pulses.BEAM==cand.BEAM) & (pulses.SAP==cand.SAP)].shape[0] * 100)
    row = [cand.id, date, vers, idL, cand.SAP, cand.BEAM, cand.N_pulses, RA, DEC, cand.DM, '{:.1f}'.format(cand.Sigma), cand_percent, '', '', 'ToProcess', '', link]    
    wks.append_row(row)  #Possible to append all together (faster)?

  #Sort spreadsheet
  wks = sh.worksheet("Classifier_v4")
  col_size = wks.col_count
  col_size_letter = chr(col_size-1 + ord('A'))
  row_size = wks.row_count
  old_cells = wks.range('A1:{col}{row}'.format(row=row_size,col=col_size_letter))
  new_cells = list(old_cells)
  old_col = wks.col_values(1)
  new_idx = [old_col.index(i) for i in sorted(old_col, reverse=True)]
  #Move repeated at the end
  format_idx = []
  temp = []
  for n in new_idx:
    if n in format_idx: temp.append(n)
    else: format_idx.append(n)
  format_idx.extend(temp)
  for i, i_new in enumerate(format_idx):
    for col in range(col_size):
      new_cells[i*col_size+col].value = old_cells[i_new*col_size+col].input_value
  wks.update_cells(new_cells)

  #Resize spreadsheet
  try:
    row = wks.col_values(1)
    row = [val for val in row if val != '']
    row = len(row)
    wks.resize(rows=row, cols=col_size)
  except ValueError: pass

  return
