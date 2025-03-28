.. _cgc-information-generation:

CAZyme Gene Cluster (CGC) Annotation
====================================

Introduction
-------------

A CAZyme Gene Cluster (CGC) refers to a group of genes co-located on the genome that collectively participate in glycan metabolism. These clusters encode various enzymes and regulatory proteins that work together to process specific carbohydrate substrates. CGCs are particularly important in microbial genomes, where they enable efficient utilization of diverse carbohydrate sources in the environment.

Generating CGC Information
---------------------------

To identify and analyze CGCs, users must first prepare annotation information by converting their GFF file into a specialized CGC-ready format. The ``gff_process`` command handles this conversion through several steps:

1. Extracts non-CAZyme sequences from protein files
2. Uses DIAMOND to annotate Transporters (TCs)
3. Employs pyHMMER to identify Transcription Factors (TFs) and Signal Transduction Proteins (STPs)
4. Combines and filters all results based on coverage and e-value thresholds
5. Integrates these annotations into the user-submitted GFF file to generate a ``cgc.gff`` file

Input GFF File Types
-----------------------

The processing workflow varies depending on the source and organism type of your GFF file:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - GFF Source
     - Processing Details
   * - NCBI Prokaryotic(``--gff_type NCBI_prok``)
     - Extracts CDS corresponding to genes to obtain protein IDs and merges this information
   * - Prodigal-generated(``--gff_type prodigal``)
     - Directly processes CDS information without additional extraction
   * - NCBI Eukaryotic(``--gff_type NCBI_euk``)
     - Extracts mRNA and CDS information for coding genes; handles non-coding genes separately (beta feature)
   * - JGI(``--gff_type JGI``)
     - Extracts protein IDs corresponding to genes from JGI-formatted annotation

Command Examples
------------------

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff EscheriaColiK12MG1655.gff --gff_type NCBI_prok

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff output_dir/uniInput.gff --gff_type prodigal

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff Xylhe1_GeneCatalog_proteins.gff --gff_type JGI

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff eukaryotic_genome.gff --gff_type NCBI_euk

.. note::

   Eukaryotic GFF processing is currently in beta. While functional, we're still validating and improving this feature.

Key Parameters
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - Directory for both input and output files, it is generated from the ``CAZyme annotation`` step
   * - ``--db_dir``
     - Directory containing database files, it is generated from the ``database`` step
   * - ``--input_gff``
     - Path to your GFF annotation file
   * - ``--gff_type``
     - Format of the GFF file: ``NCBI_prok``, ``prodigal``, ``JGI``, or ``NCBI_euk``


Output Files
-------------

The CGC finder generates several output files:

* ``diamond.out.tc`` - DIAMOND output for Transporters (TCs) via TCDB database (https://www.tcdb.org/)
* ``TF_hmm_results.tsv`` - pyHMMER output for Transcription Factors (TFs) using the TF-HMM database.
* ``STP_hmm_results.tsv`` - pyHMMER output for Signal Transduction Proteins (STPs) using the STP-HMM database.
* ``total_cgc_info.tsv`` - Comprehensive table containing all CGC information, including TCs, TFs, and STPs.

.. admonition:: Next Steps

   After generating the CGC annotation information, proceed to :doc:`CGC identification <CGC_annotation>` to identify and analyze CAZyme gene clusters in your genome.




