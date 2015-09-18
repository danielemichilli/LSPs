import subprocess
import os






def upload(folder):
  
  
  pass

  path = os.path.dirname(folder)
  subprocess.call(['sh',os.path.dirname(os.path.realpath(__file__))+'/p.sh','ls','-1'])  #MODIFICARE!!
  
  #AGGIUNGERE test per verificare che cartelle nel sito, lista nel sito, lista in sheet e files su cartesius coincidono
  ##Compare file list on website and sheet
  #site_list = subprocess.check_output(["sh", "/home/danielem/scripts/LSPs_file_list.sh"]).split('\n')[3:-1]
  #site_list = sorted(site_list)
  #sheet_list = 
  
  
  
  
  
  
  #upgrade html list file
  #create html file in the candidates folder
  #send candidates folder on the website
  #remove html file in the candidates folder
  #upgrade sheet
  
  return
