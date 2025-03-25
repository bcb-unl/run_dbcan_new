import os
import pytest
from click.testing import CliRunner
from pathlib import Path

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


# Integration test with database creation in test
def test_easy_substrate_integration(runner, tmp_path):
    """
    Integration test for easy_substrate command with database created within the test.
    """
    # Create temporary directories
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Convert to string paths for commands
    db_dir_str = str(db_dir)
    output_dir_str = str(output_dir)

    # Print test setup info
    print(f"Setting up test database in: {db_dir_str}")

    # Build database directly in test
    db_result = runner.invoke(cli, [
        'database',
        '--db_dir', db_dir_str
    ])

    # Check if database command succeeded
    if db_result.exit_code != 0:
        print(f"Database build failed: {db_result.output}")
        print(f"Exception: {db_result.exception}")
        pytest.skip("Failed to build database, skipping test")

    # Verify test files exist
    assert os.path.exists(TEST_PROTEIN), f"Test protein file not found at {TEST_PROTEIN}"
    assert os.path.exists(TEST_GFF), f"Test GFF file not found at {TEST_GFF}"

    # Print test information for debugging
    print(f"Running test with:")
    print(f"  TEST_PROTEIN: {TEST_PROTEIN}")
    print(f"  TEST_GFF: {TEST_GFF}")
    print(f"  db_dir: {db_dir_str}")
    print(f"  output_dir: {output_dir_str}")

    # Run actual command (without mocking)
    result = runner.invoke(cli, [
        'easy_substrate',
        '--mode', 'protein',
        '--input_raw_data', TEST_PROTEIN,
        '--input_gff', TEST_GFF,
        '--gff_type', 'NCBI_prok',
        '--output_dir', output_dir_str,
        '--db_dir', db_dir_str
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


