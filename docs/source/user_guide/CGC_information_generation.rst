Prepare the CGC annotation information
===========================================

A CAZyme gene cluster (CGC) refers to a group of genes co-located on the genome that are collectively involved in the metabolism of carbohydrates. These gene clusters encode enzymes and other proteins that work together to perform specific functions related to carbohydrate processing. The concept of a CAZyme gene cluster is particularly relevant in the context of microbial genomes, where such clusters often play crucial roles in the utilization of diverse carbohydrate sources.

If users want to predict CGCs, this step will be necessary to focuses on theconvert user submitted gff file into the cgc.gff generationfile. First, it extracteds the non-CAZyme sequences from faa file, and then applied runs DIAMOND to annotate TCs. Next it useds pyhmmer to annotate TFs and STPs. All three results were are combined to and filtered to keep the best hits based on the coverage and evalue (same as dbCAN domain filter rules). Annotations were are labeled added in the user submitted gff file to generate the cgc.gff.

For prokaryotic gff:

If downloaded from NCBI mode: extract the CDS corresponding to the gene to obtain the protein ID, and add it to the gene to obtain the gene information and the corresponding protein ID.

If predicted with prodigal mode: directly obtain CDS information.



For eukaryotic gff (beta function, validation is ongoing):


If downloaded from NCBI: extract the mRNA and CDS information corresponding to the gene, obtain the protein ID and add it to the gene; for non-coding genes (such as tRNA), extract the corresponding information and add it to the gene. Finally, the cgc.gff of the complete gene is obtained, which contains all the gene locations and corresponding functions.

If downloaded from JGI: extract the protein ID corresponding to the gene.


.. code-block:: shell

    run_dbcan gff_process --output_dir output_EscheriaColiK12MG1655_faa --db_dir db --input_gff EscheriaColiK12MG1655.gff --gff_type NCBI_prok

.. code-block:: shell

    run_dbcan gff_process --output_dir output_EscheriaColiK12MG1655_fna/ --db_dir db --input_gff  output_EscheriaColiK12MG1655_fna/uniInput.gff --gff_type prodigal

.. code-block:: shell

    run_dbcan gff_process --output_dir output_Xylona_heveae_TC161_faa/ --db_dir db --input_gff Xylhe1_GeneCatalog_proteins_20130827.gff --gff_type JGI

.. code-block:: shell

    run_dbcan gff_process --output_dir output_Xylhe1_faa/ --db_dir db --input_gff Xylona_heveae_TC161.gff  --gff_type NCBI_euk




