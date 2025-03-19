User Guide
==========

Update: This is the new version of run_dbCAN. We add multiple new features and improve the performance of the pipeline. The new version is more user-friendly and more efficient.
We recommend users to use the new version of run_dbCAN. If you have any questions or suggestions, please feel free to contact us.

All conda environments can be found at https://github.com/bcb-unl/run_dbcan_new/tree/master/envs

1. Add a function for downloading database files which is simpler than before.

2. Import pyrodigal (https://pyrodigal.readthedocs.io/en/stable/) instead of prodigal for input processing, besides, add function for data preprocessing, and now  run_dbCAN could support prodigal format, JGI format, and NCBI format with setting parameter.

3. Import pyHMMER (https://pyhmmer.readthedocs.io/en/stable/) instead of HMMER, which is more efficient and speeds up than HMMER.  Redesigned memory usage, now can use less memory, or high memory + high efficiency.

4. Re-organized the logic and structure of run_dbCAN. Now we split functions into each module and use “CLASS” to handle it, which is easier to update and control. Besides, use python to rewrite almost non-python codes and it's more readable.  Use config to organize all parameters.

5. Use pandas for data processing.

6. Add coverage justifications and location information in dbCAN-sub.

7. Add CAZyme justification in the final result (extra column called "Best Results).

8. Added a lot of log processing and time reporting, making it more user-friendly.

9. Re-design the CGCFinder (Now support JGI, NCBI, prodigal formats, and could directly search eukaryotes such as fungi genomes).

10. Change the blastp search to DIAMOND search in substrate prediction part, which is faster and more efficient.

11. Update steps for metagenomic data protocols.

Hint:If you want to run from raw reads from metagenome, please refer to Run from Raw Reads: Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol.
Otherwise, please refer to any following instruction. Please note that some of the precomputed results have different names from the previous version.



.. toctree::
   :maxdepth: 1
   :caption: getting_started

   getting_started/installation
   getting_started/quick_start



.. toctree::
   :maxdepth: 1
   :caption: user_guide

   user_guide/prepare_the_database
   user_guide/CAZyme_annotation
   user_guide/CGC_information_generation
   user_guide/CGC_annotation
   user_guide/predict_CGC_substrate
   user_guide/CGC_plots


.. toctree::
   :maxdepth: 1
   :caption: API

   api/index
   api/database
   api/CAZyme_annotation
   api/gff_process
   api/cgc_finder
   api/substrate_prediction
   api/cgc_circle_plot
   api/easy_workflow

.. toctree::
   :maxdepth: 1
   :caption: metagenomics_pipeline

   metagenomics_pipeline/run_from_raw_reads
   metagenomics_pipeline/run_from_raw_reads_am
   metagenomics_pipeline/run_from_raw_reads_pr
   metagenomics_pipeline/run_from_raw_reads_wk
   metagenomics_pipeline/run_from_raw_reads_em
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_co_assem
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_subsample
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_assem_free
