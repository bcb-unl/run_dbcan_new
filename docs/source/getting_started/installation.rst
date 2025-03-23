Installation
=============


We currently support configuration using anaconda and pip (docker will be published as soon as possible).

.. hint::
  We recommend using conda to install the package, as it is more stable and easier to manage dependencies.
  If you want to use pip, please make sure that all dependencies are installed correctly (have to download diamond manually because pypi doesn't support it).


.. highlights::
  - conda
  - pip
  - docker (coming soon)


.. hint::
  we didn't ask the users to git clone all repo from the github because we've uploaded it to the Pypi, and the users can install it by pip.
  Users need to prepare the environmental files, which could be downloaded from the github repo. We provide the environmental files in the `envs` folder, and could also be found directly in this link:
  prepare the conda environment (available at https://github.com/bcb-unl/run_dbcan_new/tree/master/envs)

.. code-block:: shell

  conda env create -f environment.yml
  conda activate run_dbcan


.. hint::

  Please double check the conda version, if you have a lower version of conda, please update it to the latest version.

  and make sure your python/pip path is the same as the conda path.
  If you have multiple python/pip versions, please check the version of python/pip in your conda environment.


.. code-block:: shell

  conda --version

  conda update --all

  which python

  which pip

or users can use pip to install the package (need to install `diamond` in the environment first):

.. code-block:: shell

  pip install dbcan==5.0.0
