Invoking APBS
=============

Unless invoked by :doc:`/pdb2pqr/index` or another software package (see :doc:`/other-software`), most APBS calculations are performed from the command line by running an APBS input file through the APBS program.

As mentioned in the installation and availability section, the main APBS binary is installed in ``${APBS_PREFIX}/bin`` where ``${APBS_PREFIX}`` is the top-level directory you chose for the installation. 
You can move the binary to any directory you choose.

APBS is invoked with a very simple syntax:

.. code-block:: bash
   
   apbs [options] input-file

where the list of ``[options]`` can be obtained by running APBS with the ``--help`` option.
The input file format is described in :doc:`input/index`.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   input/index.rst