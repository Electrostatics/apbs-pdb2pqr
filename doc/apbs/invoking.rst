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

The `Opal Toolkit <http://nbcr.ucsd.edu/data/docs/opal/>`_ is software produced by `NBCR <http://nbcr.ucsd.edu/wordpress2/>`_.
This toolkit allows for the computing load for processor intensive scientiﬁc applications to be shifted to a 3rd party and/or generic computing grid.
This can be tremendously advantageous in situations where a large amount of computing power is not locally available, but is required, for the task at hand.
In particular, many users have discovered that their local computational resources are insufficient for certain types of APBS calculations on large systems or at extremely high accuracy.
This client removes this resource limitation by allowing users to run on clusters at NBCR.
Recent developmental versions APBS add optional support for the off-loading of APBS calculations to an Opal service. 
Opal support has been integrated into APBS such that the end user will not be able to tell the difference between local and Opal runs of APBS: the APBS Opal client can be invoked in exactly the same way as the main APBS binary with identical output.
The APBS Opal support is in the form of a Python script ``ApbsClient.py`` and is installed by default.
The script has been tested on Python 2.5; newer/older versions of Python may or may be functional.
As mentioned above, the basic invocation is the same as the main binary.

.. code-block:: bash
   
   ApbsClient.py [options] input-file

.. toctree::
   :maxdepth: 2
   :caption: Contents

   input/index.rst