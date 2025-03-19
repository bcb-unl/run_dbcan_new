Run with CGCFinder
==================

Once cgc.gff(saved in the output_dir) is created, users can use it to predict CGCs:


.. code-block:: shell

    run_dbcan cgc_finder --output_dir output_EscheriaColiK12MG1655_faa

.. code-block:: shell

    run_dbcan cgc_finder --output_dir output_EscheriaColiK12MG1655_fna/

.. code-block:: shell

    run_dbcan cgc_finder --output_dir output_Xylona_heveae_TC161_faa/

.. code-block:: shell

    run_dbcan cgc_finder --output_dir output_Xylhe1_faa/


Noted: Now we support the following two rules to predict CGCs:

1. Forward and backward search with num of the null genes, once it finds the core/additional genes, extend to the next loop.

2. Following the AntiSMASH, we use base-pair distance (default as 15k) to search forward and backward core/additional genes, the base-pair distance is defined as the distance between two sig genes.


Users could use rule 1 (default) or rule 2,  or apply both to limit a stricter CGC.

.. code-block:: shell

    run_dbcan CGC_annotation --output_dir output_EscheriaColiK12MG1655_faa --use_null_genes true --num_null_gene 5 --use_distance true --base_pair_distance 15000 --additional_genes TC TF STP
