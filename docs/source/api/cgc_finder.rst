CGCFinder
==================

.. click:: dbcan.main:cli
   :prog: run_dbcan
   :commands: cgc_finder
   :nested: full


CGCFinder is a tool to identify carbohydrate-active enzyme (CAZyme) gene clusters (CGCs) in a genome. It uses a combination of CAZyme annotation and clustering algorithms to identify CGCs based on the presence of CAZymes and their genomic context.


Usage Examples
----------------
Identify CGCs in a genome (need run CAZyme_annotation first, or could run easy_CGC directly )
.. code-block:: bash

   run_dbcan cgc_finder --output_dir test
