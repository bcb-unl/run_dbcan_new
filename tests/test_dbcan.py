import pytest
from unittest.mock import patch
from dbcan.run_dbcan import main as run_dbcan_main
import os
from pathlib import Path

TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")

@pytest.mark.parametrize("args", [
    ['run_dbcan', 'database', '--db_dir', os.path.join(DATA_ROOT, 'db')],
    ['run_dbcan', 'easy_substrate', '--input_raw_data', os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.fna'), '--mode', 'prok', '--output_dir', os.path.join(DATA_ROOT, 'output_EscheriaColiK12MG1655_fna_sub'), '--db_dir', os.path.join(DATA_ROOT, 'db'), '--input_gff', os.path.join(DATA_ROOT, 'gff'), '--input_gff_format', 'prodigal'],
    ['run_dbcan', 'easy_substrate', '--input_raw_data', os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.faa'), '--mode', 'protein', '--output_dir', os.path.join(DATA_ROOT, 'output_EscheriaColiK12MG1655_faa_sub'), '--db_dir', os.path.join(DATA_ROOT, 'db'), '--input_format', 'NCBI', '--input_gff', os.path.join(DATA_ROOT, 'EscheriaColiK12MG1655.gff'), '--input_gff_format', 'NCBI_prok'],
    ['run_dbcan', 'easy_substrate', '--input_raw_data', os.path.join(DATA_ROOT, 'Xylona_heveae_TC161.faa'), '--mode', 'protein', '--output_dir', os.path.join(DATA_ROOT, 'output_Xylona_heveae_TC161_faa_sub'), '--db_dir', os.path.join(DATA_ROOT, 'db'), '--input_format', 'NCBI', '--input_gff', os.path.join(DATA_ROOT, 'Xylona_heveae_TC161.gff'), '--input_gff_format', 'NCBI_euk'],
    ['run_dbcan', 'easy_substrate', '--input_raw_data', os.path.join(DATA_ROOT, 'Xylhe1_GeneCatalog_proteins_20130827.aa.fasta'), '--mode', 'protein', '--output_dir', os.path.join(DATA_ROOT, 'output_Xylhe1_faa_sub'), '--db_dir', os.path.join(DATA_ROOT, 'db'), '--input_format', 'JGI', '--input_gff', os.path.join(DATA_ROOT, 'Xylhe1_GeneCatalog_proteins_20130827.gff'), '--input_gff_format', 'JGI']
])
def test_run_dbcan_commands(args):
    with patch('sys.argv', args):
        run_dbcan_main()
