Substrate prediction with CGCs
==================================

Since 2023, we included a function to predict substrates for CGCs. It is based on two methods, which have been described in our dbCAN3 paper. 
Note: Change BLASTP into DIAMOND BLASTP in the substrate prediction part, which is faster and more efficient.

.. code-block:: shell

    run_dbcan CGC_substrate_prediction --input output_EscheriaColiK12MG1655_fna --db_dir db --output_dir output_EscheriaColiK12MG1655_fna

.. code-block:: shell

    run_dbcan CGC_substrate_prediction --input output_EscheriaColiK12MG1655_faa --db_dir db --output_dir output_EscheriaColiK12MG1655_faa

.. code-block:: shell

    run_dbcan CGC_substrate_prediction --input output_Xylhe1_faa --db_dir db --output_dir output_Xylhe1_faa

.. code-block:: shell

    run_dbcan CGC_substrate_prediction --input output_Xylona_heveae_TC161_faa --db_dir db --output_dir output_Xylona_heveae_TC161_faa
