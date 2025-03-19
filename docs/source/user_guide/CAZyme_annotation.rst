CAZyme Annotation
=========================

The purpose of this step is to annotate CAZymes based on the input raw data (prokaryotic genome, metagenomics contigs or protein sequences). The input data can be in FASTA format, and the output will be a table with CAZyme annotations. The table will include information such as the CAZyme family, the number of genes, and the location of the genes. The output will also include a summary of the CAZyme annotations.

Could use "--methods"  to choose different combination of tools from "diamond", "hmm", and "dbCANsub".

.. code-block:: shell

    run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.fna	 --out_dir output_EscheriaColiK12MG1655_fna --db_dir db --mode prok

.. code-block:: shell

    run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.faa   --out_dir output_EscheriaColiK12MG1655_faa --db_dir db --mode protein

.. code-block:: shell

    run_dbcan CAZyme_annotation --input_raw_data Xylona_heveae_TC161.faa --out_dir output_Xylona_heveae_TC161_faa --db_dir db --mode protein

.. code-block:: shell

    run_dbcan CAZyme_annotation --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --out_dir output_Xylhe1_faa --db_dir db --mode protein

