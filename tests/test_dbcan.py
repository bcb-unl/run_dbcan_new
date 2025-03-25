import os
import pytest
from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch, MagicMock

from dbcan.main import cli

# Test root directory and data directory
TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")

# Test data file paths
TEST_PROTEIN = os.path.join(DATA_ROOT, "EscheriaColiK12MG1655.faa")
TEST_NUCLEOTIDE = os.path.join(DATA_ROOT, "EscheriaColiK12MG1655.fna")
TEST_GFF = os.path.join(DATA_ROOT, "EscheriaColiK12MG1655.gff")


@pytest.fixture
def runner():
    """Return a Click CLI test runner"""
    return CliRunner()


# Mock database creation
@pytest.fixture
def mock_db_dir(tmp_path):
    """Create a mock database directory with necessary files"""
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    # Create mock database files
    (db_dir / "CAZy.dmnd").touch()
    (db_dir / "dbCAN.hmm").touch()
    (db_dir / "dbCAN-sub.hmm").touch()

    return str(db_dir)


# Integration test with mocked database creation
@patch('dbcan.core.run_dbCAN_input_process', return_value=None)
@patch('dbcan.core.run_dbCAN_CAZyme_annotation', return_value=None)
@patch('dbcan.core.run_dbCAN_CGCFinder_preprocess', return_value=None)
@patch('dbcan.core.run_dbCAN_CGCFinder', return_value=None)
@patch('dbcan.core.run_dbCAN_CGCFinder_substrate', return_value=None)
@patch('dbcan.core.run_dbcan_syn_plot', return_value=None)
def test_easy_substrate_cmd(mock_plot, mock_substrate, mock_cgc, mock_preprocess, mock_annotation, mock_input, runner, mock_db_dir, tmp_path):
    """
    Integration test for easy_substrate command with mocked database creation.
    """
    # Create temporary output directory
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    output_dir_str = str(output_dir)

    # Verify test files exist
    assert os.path.exists(TEST_PROTEIN), f"Test protein file not found at {TEST_PROTEIN}"
    assert os.path.exists(TEST_GFF), f"Test GFF file not found at {TEST_GFF}"

    # Print test information for debugging
    print(f"Running test with:")
    print(f"  TEST_PROTEIN: {TEST_PROTEIN}")
    print(f"  TEST_GFF: {TEST_GFF}")
    print(f"  db_dir: {mock_db_dir}")
    print(f"  output_dir: {output_dir_str}")

    # Run actual command (with mocking)
    result = runner.invoke(cli, [
        'easy_substrate',
        '--mode', 'protein',
        '--input_raw_data', TEST_PROTEIN,
        '--input_gff', TEST_GFF,
        '--gff_type', 'NCBI_prok',
        '--output_dir', output_dir_str,
        '--db_dir', mock_db_dir,
        '--threads', '4'
    ])

    # Print output if there was an error
    if result.exit_code != 0:
        print(f"Command failed with exit code {result.exit_code}")
        print(f"Output: {result.output}")
        print(f"Exception: {result.exception}")

    assert result.exit_code == 0, f"Command failed: {result.output}"

    # Verify key output files were created
    assert os.path.exists(os.path.join(output_dir_str, "overview.tsv")), "overview.tsv not found"
    assert os.path.exists(os.path.join(output_dir_str, "cgc_standard_out.tsv")), "cgc_standard_out.tsv not found"
    assert os.path.exists(os.path.join(output_dir_str, "substrate_prediction.tsv")), "substrate_prediction.tsv not found"


