Installing APBS
===============


===================
Binary installation
===================

------------------
Windows
------------------

If you are using APBS on a Windows system, you may not want to install APBS in a directory with spaces in the path name (e.g., ``C:\Program Files\``) as this can cause problems with some visualization programs.
The windows release is for Windows 7 and is 64-bit only.
It installs to ``C:\APBS`` by default and provides a batch file that creates a shell with the correct path.

--------
Mac OS X
--------

The OS X build is built for Mavericks, available as a :file:`.dmg` (disk image) file with an APBS Application Bundle.
Just open the :file:`.dmg` file and drag the app to the :file:`Applications` folder.
Running the App Bundle opens a Terminal window with the appropriate path to the APBS binary.

-----
Linux
-----

APBS binaries are provided in compressed tar format (:file:`*.tgz`) for Linux. 
This tarball can be opened by

.. code-block:: bash

   gzip -dc apbs-{#}.{#}.{#}-{XYZ}.tgz | tar xvf -

where :file:`apbs-{#}.{#}.{#}-{XYZ}.tgz` is the downloaded tarball with ``XYZ`` as the particular architecture of the binary you downloaded and ``#.#.#`` as the version number.
This will expand into a directory called :file:`apbs-{#}.{#}.{#}-{XYZ}`.
The contents of this directory can be placed anywhere on your system that you prefer (and have access to).



===================
Source installation
===================

We recommend that most users compile APBS from our official releases; see :doc:`/downloads` for downloading instructions.
Adventurous users may want to try to compile the developmental versions on `GitHub <https://github.com/Electrostatics/apbs-pdb2pqr>`_.
In either case, detailed compile/build instructions are provided with the source code.

==================
What's in the box?
==================

bin
  contains the main APBS executable
share/apbs
  contains additional APBS-related files
doc
  the APBS programmer guide
examples
  APBS examples
tests
  the APBS test suite
tools
  useful programs to help process APBS input and output
include
  header files for building software that calls APBS
lib
  libraries for building software that calls APBS
