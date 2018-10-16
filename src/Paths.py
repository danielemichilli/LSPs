import os

OBS_FOLDER = os.environ["LOTAAS_OUTDIR"]
RAW_FOLDER = os.environ["LOTAAS_RAW1"]
PL_FOLDER = os.environ["LOTAAS_PL"]
WRK_FOLDER = '/dev/shm/wrk'
TMP_FOLDER = '/dev/shm/tmp'
RAW_TEMP = os.environ["LOTAAS_RAW2"]

#MODEL_FILE = os.path.join(os.environ["LOTAAS_PL"], "LSP/sp_ML.model")
MODEL_FILE = "/lustre4/0/lotaas2/software4/LOTAAS-Scripts_v2/LSP/sp_ML.model"