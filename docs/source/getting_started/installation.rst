Installation
=============


We currently support configuration using anaconda and pip (docker will be published as soon as possible).

.. hint::
  we didn't ask the users to git clone all repo from the github because we've uploaded it to the Pypi, and the users can install it by pip.
  Users need to prepare the environmental files, which could be downloaded from the github repo. We provide the environmental files in the `envs` folder, and could also be found directly in this link:
  prepare the conda environment (available at https://github.com/bcb-unl/run_dbcan_new/tree/master/envs)

.. code-block:: shell

  conda env create -f environment.yml
  conda activate run_dbcan

or users can use pip to install the package (need to install `diamond` in the environment first):
.. code-block:: shell

  pip install dbcan==5.0.0
