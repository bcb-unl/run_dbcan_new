Substrate prediction with CGCs
==================================

Since 2023, we included a function to predict substrates for CGCs. It is based on two methods, which have been described in our dbCAN3 paper.
Note: Change BLASTP into DIAMOND BLASTP in the substrate prediction part, which is faster and more efficient.

.. code-block:: shell

    run_dbcan substrate_prediction --output_dir output_EscheriaColiK12MG1655_fna --db_dir db

.. code-block:: shell

    run_dbcan substrate_prediction --output_dir output_EscheriaColiK12MG1655_faa --db_dir db

.. code-block:: shell

    run_dbcan substrate_prediction --output_dir output_Xylhe1_faa --db_dir db

.. code-block:: shell

    run_dbcan substrate_prediction --output_dir output_Xylona_heveae_TC161_faa --db_dir db
