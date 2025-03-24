import os
import shutil
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


# Fixtures
@pytest.fixture
def runner():
    """Return a Click CLI test runner"""
    return CliRunner()


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return str(output_dir)


@pytest.fixture
def temp_db_dir(tmp_path):
    """Create a temporary database directory"""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    return str(db_dir)


# Test easy_substrate command (encompasses all other functionality)
@patch('dbcan.core.run_dbCAN_input_process')
@patch('dbcan.core.run_dbCAN_CAZyme_annotation')
@patch('dbcan.core.run_dbCAN_CGCFinder_preprocess')
@patch('dbcan.core.run_dbCAN_CGCFinder')
@patch('dbcan.core.run_dbCAN_CGCFinder_substrate')
@patch('dbcan.core.run_dbcan_syn_plot')
def test_easy_substrate_cmd(mock_plot, mock_substrate, mock_cgc, mock_preprocess, mock_annotation, mock_input,
                           runner, temp_output_dir, temp_db_dir):
    """Test if easy_substrate command correctly calls underlying functions"""
    # Verify test file exists
    assert os.path.exists(TEST_PROTEIN), f"Test protein file not found at {TEST_PROTEIN}"

    result = runner.invoke(cli, [
        'easy_substrate',
        '--mode', 'protein',
        '--input_raw_data', TEST_PROTEIN,
        '--input_gff', TEST_GFF,
        '--gff_type', 'NCBI_prok',
        '--output_dir', temp_output_dir,
        '--db_dir', temp_db_dir,
        '--threads', '8',
    ])

    assert result.exit_code == 0
    mock_input.assert_called_once()
    mock_annotation.assert_called_once()
    mock_preprocess.assert_called_once()
    mock_cgc.assert_called_once()
    mock_substrate.assert_called_once()
    mock_plot.assert_called_once()

    # Validate configs passed to underlying functions
    input_config = mock_input.call_args[0][0]
    assert input_config.input_file == TEST_PROTEIN
    assert input_config.out_dir == temp_output_dir

    annotation_config = mock_annotation.call_args[0][0]
    assert annotation_config.db_dir == temp_db_dir
    assert annotation_config.methods == 'all'
    assert annotation_config.threads == 4


# Integration test for GitHub workflow
@pytest.mark.integration
def test_easy_substrate_integration(runner):
    """
    Integration test for easy_substrate command using database created in GitHub workflow.
    This test will run the actual command without mocking.
    """
    # Determine the correct paths for test execution
    # In CI environment, use the current directory
    current_dir = os.getcwd()
    db_dir = os.path.join(current_dir, "db")

    # Verify database exists
    if not os.path.exists(db_dir):
        pytest.skip(f"Database directory not found at {db_dir}")

    # Create output directory
    output_dir = os.path.join(current_dir, "test_output")
    os.makedirs(output_dir, exist_ok=True)

    # Verify test files exist
    assert os.path.exists(TEST_PROTEIN), f"Test protein file not found at {TEST_PROTEIN}"
    assert os.path.exists(TEST_GFF), f"Test GFF file not found at {TEST_GFF}"

    # Print test information for debugging
    print(f"Running test with:")
    print(f"  TEST_PROTEIN: {TEST_PROTEIN}")
    print(f"  TEST_GFF: {TEST_GFF}")
    print(f"  db_dir: {db_dir}")
    print(f"  output_dir: {output_dir}")

    # Run actual command (without mocking)
    result = runner.invoke(cli, [
        'easy_substrate',
        '--mode', 'protein',
        '--input_raw_data', TEST_PROTEIN,
        '--input_gff', TEST_GFF,
        '--gff_type', 'NCBI_prok',
        '--output_dir', output_dir,
        '--db_dir', db_dir,
        '--threads', '8'
    ])

    # Print output if there was an error
    if result.exit_code != 0:
        print(f"Command failed with exit code {result.exit_code}")
        print(f"Output: {result.output}")
        print(f"Exception: {result.exception}")

    assert result.exit_code == 0, f"Command failed: {result.output}"

    # Verify key output files were created
    assert os.path.exists(os.path.join(output_dir, "overview.tsv")), "overview.tsv not found"
    assert os.path.exists(os.path.join(output_dir, "cgc_standard_out.tsv")), "cgc_standard_out.tsv not found"
    assert os.path.exists(os.path.join(output_dir, "substrate_prediction.tsv")), "substrate_prediction.tsv not found"


