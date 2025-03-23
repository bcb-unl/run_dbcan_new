Gff procress
==================

.. click:: dbcan.main:cli
   :prog: run_dbcan
   :commands: gff_process
   :nested: full


Example Usage
----------------
.. hint::
   Need to run CAZyme annotation first. This is the subprogress of easy_CGC/easy_substrate.

.. code-block:: bash

   run_dbcan gff_process --db_dir db --input_gff test.gff --output_dir test --gff_type JGI/NCBI_prok/NCBI_euk/prodigal (choose your format)
