Substrate Prediction
====================

.. click:: dbcan.main:cli
   :prog: run_dbcan
   :commands: substrate_prediction --help
   :nested: full


example usage
----------------

.. hint::
   This is the subprocess of easy_substrate.
   Need to run CAZyme annotation and cgcfinder first.
   Or could use easy_substrate directly.


.. code:: bash
   run_dbcan substrate_prediction --output_dir test --db_dir db
