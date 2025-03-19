import os
from pathlib import Path

TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")

INPUT_PROTEIN_NAME="uniInput.faa"
INPUT_PRODIGAL_GFF_NAME="uniInput.gff"


CAZY_COLUMN_NAMES = [
    'Gene ID',
    'CAZy ID',
    '% Identical',
    'Length',
    'Mismatches',
    'Gap Open',
    'Gene Start', 
    'Gene End', 
    'CAZy Start', 
    'CAZy End', 
    'E Value', 
    'Bit Score'
]

TCDB_COLUMN_NAMES = [
    'TCDB ID', 
    'TCDB Length', 
    'Target ID', 
    'Target Length', 
    'EVALUE', 
    'TCDB START',
    'TCDB END', 
    'QSTART', 
    'QEND', 
    'COVERAGE'
]

HMMER_COLUMN_NAMES = [
    'HMM Name', 
    'HMM Length', 
    'Target Name', 
    'Target Length', 
    'i-Evalue',
    'HMM From', 
    'HMM To', 
    'Target From', 
    'Target To', 
    'Coverage',
    'HMM File Name'
]

dbCAN_sub_COLUMN_NAMES = [
    'Subfam Name', 
    'Subfam Composition', 
    'Subfam EC', 
    'Substrate', 
    'HMM Length',
    'Target Name', 
    'Target Length', 
    'i-Evalue', 
    'HMM From', 
    'HMM To', 
    'Target From',
    'Target To', 
    'Coverage', 
    'HMM File Name'
]

OVERLAP_RATIO_THRESHOLD = 0.5

CAZY_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/CAZy.dmnd"
HMMER_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN.hmm"
DBCAN_SUB_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN_sub.hmm"
DBCAN_SUB_MAP_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/fam-substrate-mapping.tsv"

TCDB_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/tcdb.dmnd"
TF_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/TF.hmm"
STP_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/STP.hmm"

PUL_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/PUL.dmnd"
PUL_MAP_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN-PUL.xlsx"
PUL_ALL_URL ="https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN-PUL.tar.gz"



FILE_PATHS = {
    'diamond': 'diamond_results.tsv',
    'dbcan_sub': 'dbCAN-sub.substrate.tsv',
    'dbcan_hmm': 'dbCAN_hmm_results.tsv'
}

COLUMN_NAMES = {
    'diamond': ['Gene ID', 'CAZy ID'],
    'dbcan_sub': ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue'],
    'dbcan_hmm': ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']
}

OVERVIEW_COLUMNS = ['Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results']

GFF_COLUMNS = ['Contig ID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

