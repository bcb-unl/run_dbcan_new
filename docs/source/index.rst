User Guide
==========

Update: What's New in run_dbCAN
================================

The new version of **run_dbCAN** introduces multiple new features and significant performance improvements, making the pipeline more user-friendly and efficient. We highly recommend users to upgrade to this version. If you have any questions or suggestions, please feel free to [contact us](mailto:support@dbcan.org).

All conda environments can be found at the following link:
`run_dbCAN Conda Environments <https://github.com/bcb-unl/run_dbcan_new/tree/master/envs>`_

Key Features and Improvements
-----------------------------

1. **Simplified Database Downloading**
   Added a new function for downloading database files, making the process simpler than before.

2. **Enhanced Input Processing**
   - Replaced `prodigal` with `pyrodigal` ([documentation](https://pyrodigal.readthedocs.io/en/stable/)) for input processing.
   - Added support for multiple input formats, including Prodigal, JGI, and NCBI formats, with configurable parameters.

3. **Improved HMMER Performance**
   - Replaced `HMMER` with `pyHMMER` ([documentation](https://pyhmmer.readthedocs.io/en/stable/)), which is faster and more efficient.
   - Redesigned memory usage to support both low-memory and high-efficiency modes.

4. **Modular Code Structure**
   - Reorganized the logic and structure of `run_dbCAN` by splitting functions into modules and using Python classes for better maintainability.
   - Rewrote non-Python code in Python for improved readability.
   - Centralized parameter management using configuration files.

5. **Data Processing with Pandas**
   Leveraged the power of `pandas` for efficient data processing.

6. **Enhanced dbCAN-sub Features**
   - Added coverage justifications and location information.
   - Included CAZyme justification in the final results with an extra column called "Best Results."

7. **Improved Logging and Time Reporting**
   Added extensive logging and time reporting to make the pipeline more user-friendly.

8. **Redesigned CGCFinder**
   - Now supports JGI, NCBI, and Prodigal formats.
   - Directly searches eukaryotic genomes, including fungi.

9. **Faster Substrate Prediction**
   Replaced `blastp` with `DIAMOND` for substrate prediction, significantly improving speed and efficiency.

10. **Updated Metagenomic Protocols**
    Improved steps for metagenomic data processing.

**Hint:**
If you want to run the pipeline from raw metagenomic reads, please refer to the following guide:
*Run from Raw Reads: Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol.*

Otherwise, refer to the instructions below. Please note that some precomputed results may have different names compared to the previous version.

.. note::
   For detailed instructions, refer to the respective sections in the documentation.

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
