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


@pytest.fixture
def actual_db_dir(tmp_path, runner):
    """Create an actual database directory by running the database command"""
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    # Run the database command to create the database
    result = runner.invoke(cli, [
        'database',
        '--db_dir', str(db_dir)
    ])

    if result.exit_code != 0:
        print(f"Database command failed with exit code {result.exit_code}")
        print(f"Output: {result.output}")
        print(f"Exception: {result.exception}")
        pytest.skip("Failed to create database, skipping test")

    return str(db_dir)


class TestEasySubstrate:
    @patch('dbcan.core.run_dbCAN_input_process', return_value=None)
    @patch('dbcan.core.run_dbCAN_CAZyme_annotation', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder_preprocess', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder_substrate', return_value=None)
    @patch('dbcan.core.run_dbcan_syn_plot', return_value=None)
    def test_easy_substrate_protein(self, mock_plot, mock_substrate, mock_cgc, mock_preprocess, mock_annotation, mock_input, runner, actual_db_dir, tmp_path):
        """
        Integration test for easy_substrate command with protein input.
        """
        # Create temporary output directory
        output_dir = tmp_path / "output_protein"
        output_dir.mkdir()
        output_dir_str = str(output_dir)

        # Verify test files exist
        assert os.path.exists(TEST_PROTEIN), f"Test protein file not found at {TEST_PROTEIN}"
        assert os.path.exists(TEST_GFF), f"Test GFF file not found at {TEST_GFF}"

        # Print test information for debugging
        print(f"Running test with:")
        print(f"  TEST_PROTEIN: {TEST_PROTEIN}")
        print(f"  TEST_GFF: {TEST_GFF}")
        print(f"  db_dir: {actual_db_dir}")
        print(f"  output_dir: {output_dir_str}")

        # Run actual command (with mocking)
        result = runner.invoke(cli, [
            'easy_substrate',
            '--mode', 'protein',
            '--input_raw_data', TEST_PROTEIN,
            '--input_gff', TEST_GFF,
            '--gff_type', 'NCBI_prok',
            '--output_dir', output_dir_str,
            '--db_dir', actual_db_dir,
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

    @patch('dbcan.core.run_dbCAN_input_process', return_value=None)
    @patch('dbcan.core.run_dbCAN_CAZyme_annotation', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder_preprocess', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder', return_value=None)
    @patch('dbcan.core.run_dbCAN_CGCFinder_substrate', return_value=None)
    @patch('dbcan.core.run_dbcan_syn_plot', return_value=None)
    def test_easy_substrate_nucleotide(self, mock_plot, mock_substrate, mock_cgc, mock_preprocess, mock_annotation, mock_input, runner, actual_db_dir, tmp_path):
        """
        Integration test for easy_substrate command with nucleotide input.
        """
        # Create temporary output directory
        output_dir = tmp_path / "output_nucleotide"
        output_dir.mkdir()
        output_dir_str = str(output_dir)

        # Verify test files exist
        assert os.path.exists(TEST_NUCLEOTIDE), f"Test nucleotide file not found at {TEST_NUCLEOTIDE}"

        # Print test information for debugging
        print(f"Running test with:")
        print(f"  TEST_NUCLEOTIDE: {TEST_NUCLEOTIDE}")
        print(f"  db_dir: {actual_db_dir}")
        print(f"  output_dir: {output_dir_str}")

        # Run actual command (with mocking)
        result = runner.invoke(cli, [
            'easy_substrate',
            '--mode', 'prok',
            '--input_raw_data', TEST_NUCLEOTIDE,
            '--input_gff', os.path.join(output_dir_str, "uniInput.gff"),
            '--gff_type', 'prodigal',
            '--output_dir', output_dir_str,
            '--db_dir', actual_db_dir,
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


