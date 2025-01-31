Database Installation Command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
prepare the conda environment:

.. code-block:: shell

    conda env create -f environment.yml (available at https://github.com/Xinpeng021001/run_dbCAN_new/blob/master/environment.yml)
    conda activate dbCAN-test-env
    pip install run-dbcan-new


1. We provide command in the run_dbcan script:

.. code-block:: shell

    run_dbcan database --db_dir db

2. Users could also download all  database files from the dbCAN2 website (http://bcb.unl.edu/dbCAN2/download/Databases/), and then put them into the db directory.

.. code-block:: shell

    wget -q https://bcb.unl.edu/dbCAN2/download/test/dbCAN_db.tar.gz -O db.tar.gz
    tar -zxvf db.tar.gz
    rm db.tar.gz
