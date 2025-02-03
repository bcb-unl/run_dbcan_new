Process the input fasta file
=================================

The purpose of this step is to process the input FASTA file:
For the prokaryote nucleotide sequences or metagenomic sequences, pyrodigal is applied to predict the protein-coding genes.
For the protein sequence, the main idea is to reformat ID structure to avoid potential mapping issues.

We provide four types of input data to test (all files are saved in the `example folder`_): 

1. NCBI prokaryotic genome data (Ecoli genome fastafna) 

.. code-block:: shell

    run_dbcan input_process --input_raw_data EscheriaColiK12MG1655.fna --mode prok --output_dir output_EscheriaColiK12MG1655_fna --db_dir db

2. NCBI prokaryotic protein data (Ecoli faa) 

.. code-block:: shell

    run_dbcan input_process --input_raw_data EscheriaColiK12MG1655.faa --mode protein --output_dir output_EscheriaColiK12MG1655_faa --db_dir db

3. NCBI eukaryotic protein data (fungi) 

.. code-block:: shell

    run_dbcan input_process --input_raw_data Xylona_heveae_TC161.faa --mode protein --output_dir output_Xylona_heveae_TC161_faa --db_dir db

4. JGI eukaryotic protein data (fungi) 

.. code-block:: shell

    run_dbcan input_process --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --mode protein --output_dir output_Xylhe1_faa --db_dir db

For the eukaryotic genome sequence, we suggest users using other tools to predict the protein sequence such as:

Funannotate: https://funannotate.readthedocs.io/en/latest/

BRAKER:      https://github.com/Gaius-Augustus/BRAKER

EviAnn:      https://github.com/alekseyzimin/EviAnn_release

.. _example folder: https://bcb.unl.edu/dbCAN2/download/test
