import os

OBS_FOLDER = os.environ["LOTAAS_OUTDIR"]
RAW_FOLDER = os.environ["LOTAAS_RAW1"]
SITE_CERT = os.path.join(os.environ["LOTAAS_PL"], "LSP/.site_cert.json")
WRK_FOLDER = '/dev/shm/wrk/{}'
TMP_FOLDER = '/dev/shm/tmp/{}'
MASK_FILE = os.path.join(os.environ["LOTAAS_RAW2"], '{idL}_red/SAP{sap}/BEAM{beam}/{idL}_SAP{sap}_BEAM{beam}_rfifind.mask')

MODEL_FILE = '/lustre4/0/lotaas2/software/LSP/sp_ML.model'
CLASSIFIER = '/lustre4/0/lotaas2/software/LOTAAS-Scripts/scores_robLyon/PulsarProcessingScripts-master/ML.jar'