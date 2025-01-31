Quick Start
===========

This section provides a quick guide to running the run_dbcan tool suite with example data and explains the output files generated.

For the updated run_dbcan, we provide two types of approach for users:

1.Automated analysis that can be done with one line of command.

2.One by one command that allows checking problems by breaking down the steps or making autonomous changes to some of the results.

Here we show all the steps in one line of command.

We provide multiple example data sets for users to test the tool suite. The example data sets are available in the `example_data` directory.


We currently support configuration using anaconda and pip:

.. code-block:: shell
  
  conda env create -f environment.yml (available at https://github.com/Xinpeng021001/run_dbCAN_new/blob/master/environment.yml)
  conda activate dbCAN-test-env
  pip install run-dbcan-new


1. Running Example Data for CAZyme Annotation
-----------------------------------------------

To run the dbCAN tool suite on the `Escherichia coli Strain MG1655`_ example data, use the following command. The input file `EscheriaColiK12MG1655.fna` represents the FASTA format complete genome DNA sequence.

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data EscheriaColiK12MG1655.fna --mode prok --output_dir output_EscheriaColiK12MG1655_fna --db_dir db

.. _Escherichia coli Strain MG1655: https://www.ncbi.nlm.nih.gov/nuccore/U00096.2

For the protein sequence input, use the following command (Please note that the input format is needed for the protein sequence only. `NCBI` represents the fasta ID format like NCBI ">WP_000002088.1", and the `JGI` mode represents the fasta ID format like JGI ">jgi|Xylhe1|242238|").:

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data EscheriaColiK12MG1655.faa --mode protein --output_dir output_EscheriaColiK12MG1655_faa --db_dir db --input_format NCBI

We also provide eukaryotes example data sets. For example, to run the dbCAN tool suite on the `Xylona heveae TC161`_ example data, use the following command:

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data Xylona_heveae_TC161.faa --mode protein --output_dir output_Xylona_heveae_TC161_faa --db_dir db  --input_format NCBI
.. _Xylona heveae TC161: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001619985.1/

And JGI dataset `Xylhe1`_:

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --mode protein --output_dir output_Xylhe1_faa --db_dir db  --input_format JGI

.. _Xylhe1: https://mycocosm.jgi.doe.gov/Xylhe1/Xylhe1.home.html

2. Understanding the Output
---------------------------

After running the tool, several output files are generated in the output folder, each with specific information:

**uniInput.faa**
  The unified input file for subsequent tools, created by Prodigal if a nucleotide sequence is used, or provided by the user as protein sequence.

**dbCAN-sub.substrate.tsv**
  Output from the pyHMMER using dbCAN_sub-HMM.

**diamond_results.tsv**
  Results from the Diamond BLAST using CAZy.faa.

**dbCAN_hmm_results.tsv**
  Output from the pyHMMER using dbCAN-HMM..

**overview.tsv**
  Summarizes CAZyme predictions across tools. We recommend results using at least two tools (Shown as the "Recommend Results").



3. Running Example Data for CGC Annotation
-----------------------------------------------
.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data EscheriaColiK12MG1655.fna --mode prok --output_dir output_EscheriaColiK12MG1655_fna_CGC --db_dir db --input_gff gff --input_gff_format prodigal


.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data EscheriaColiK12MG1655.faa --mode protein --output_dir output_EscheriaColiK12MG1655_faa_CGC --db_dir db --input_format NCBI --input_gff EscheriaColiK12MG1655.gff --input_gff_format NCBI_prok


.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data Xylona_heveae_TC161.faa --mode protein --output_dir output_Xylona_heveae_TC161_faa_CGC --db_dir db  --input_format NCBI --input_gff Xylona_heveae_TC161.gff --input_gff_format NCBI_euk


.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --mode protein --output_dir output_Xylhe1_faa_CGC --db_dir db  --input_format JGI --input_gff Xylhe1_GeneCatalog_proteins_20130827.gff --input_gff_format JGI




4. Understanding the Output
---------------------------

including the output files from the previous step, and new outputs:

**non_CAZyme.faa**
  The non-CAZyme protein sequences extracted from uniInput.faa, which is based on the overview results.

**TC_results.tsv**
  Results from the Diamond BLAST using TCDB to annotate transporter protein.

**TF_results.tsv**
  Results from the pyHMMER using TF-HMM to annotate transcription factor protein.

**STP_results.tsv**
  Results from the pyHMMER using STP-HMM to annotate signal transduction protein.

**total_cgc_info.tsv**
  The total annotation of all signature proteins combing TC, TF, and STP. Using the same overlap method to filter as CAZyme annotation.

**cgc.gff**
  The input file of CGCFinder in gff format. This is generated by the tool suite based on the input_gff file and "total_cgc_info.tsv".

**cgc_standard_out.tsv**
  The standard output of CGCFinder.


5. Running Example Data for Substrate Prediction
----------------------------------------------------

.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data EscheriaColiK12MG1655.fna --mode prok --output_dir output_EscheriaColiK12MG1655_fna_sub --db_dir db --input_gff gff --input_gff_format prodigal


.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data EscheriaColiK12MG1655.faa --mode protein --output_dir output_EscheriaColiK12MG1655_faa_sub --db_dir db --input_format NCBI --input_gff EscheriaColiK12MG1655.gff --input_gff_format NCBI_prok


.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data Xylona_heveae_TC161.faa --mode protein --output_dir output_Xylona_heveae_TC161_faa_sub --db_dir db  --input_format NCBI --input_gff Xylona_heveae_TC161.gff --input_gff_format NCBI_euk


.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --mode protein --output_dir output_Xylhe1_faa_sub --db_dir db  --input_format JGI --input_gff Xylhe1_GeneCatalog_proteins_20130827.gff --input_gff_format JGI





6. Understanding the Output
---------------------------
including the output files from the previous step, and new outputs:

**substrate.out**
  The final output of substrate prediction, which includes the substrate prediction results of each CAZyme gene cluster.

**PUL_blast.out**
  The DIAMOND blastp results of CGCs against dbCAN-PULs.

**synteny_pdf/**
  The synteny plot folder including predicted results. The plot shows the gene cluster mapping between PULs and CGCs.