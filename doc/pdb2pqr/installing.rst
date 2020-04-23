Installing PDB2PQR
==================

.. note::

   Please `register <http://eepurl.com/by4eQr>`_ before using PDB2PQR!

Most users will use PDB2PQR through `the web server <http://apbs.poissonboltzmann.org>`_ (after `registering <http://eepurl.com/by4eQr>`_, of course).
However, it is also possible to install local versions of PDB2PQR.
These local installations give a command line version of the PDB2PQR software that can be customized through a variety of extensions and used as a local web server (if compiled from source).

===================
Binary installation
===================

Everything needed to run PDB2PQR is included in the binary builds (including Python interpreter, dependencies, etc).
Just unpack and use.
PDB2PQR can be run with the “pdb2pqr” executable found in the base folder of the uncompressed archive.

Binary tarballs have ``bin`` in the file name and are named by target platform.

Binary builds *do not* provide local web server functionality.
PDB2PQR must be compiled from source to provide a local web server.

Binaries are provided for the following platforms:

* OSX binaries require OSX 10.6 or newer. The OSX binary is 64-bit.

* Linux binaries require CentOS 6 or newer and have been tested on Ubuntu 12.04 LTS and Linux Mint 13. If you are running 64-bit Linux use the 64-bit libraries.

* Windows binaries are 64 bit and were built and tested on Windows 7 64-bit but should work on Windows XP, Vista, and 8 on 64-bit systems.

===================
Source installation
===================

As the bulk of the PDB2PQR code is written Python, the PDB2PQR code itself is (mostly) architecture- and compiler-independent.
PDB2PQR has been tested using Python versions 2.6-2.7; PDB2PQR will not work with older versions of Python.
Users who simply want to use the PDB2PQR without ligand parameterization support can unarchive the source code, change to the top-level source code directory, and run:

.. code-block:: bash

  $ python scons/scons.py BUILD_PDB2PKA=False 
  $ python scons/scons.py install

If NumPy is unavailable, PDB2PQR will be built without ligand or PDB2PKA [#PDB2PKA]_ support.

.. [#PDB2PKA] PDB2PKA is the PDB2PQR library that includes both ligand parameterization and Poisson-Boltzmann-based pKa calculation routines. This code is written in C++ and Python. This portion of the code also requires the Python NumPy package.

----------------------------
Configuring the installation
----------------------------

Compilation and installation can be configured by editing the :file:`build_config.py` file.
This is the preferred way to configure the program. 
Instructions and examples for each setting are included in the file.

Configuration command-line parameters can also be used and will override any settings in :file:`build_config.py`:

``PREFIX=<DIR>``
  Set install directory. Default is :file:`~/pdb2pqr`

``URL=<URL>``
  Set url for a local webserver.  Default if ``http://<COMPUTER NAME>/pdb2pqr/``

``APBS=<APBS_BINARY>``
  Location of APBS binary.

``OPAL=<OPAL_URL>``
  Set URL for the Opal web service.

``APBS_OPAL=<APBS_OPAL_URL>``
  Set URL for the APBS Opal web service.

``MAX_ATOMS=<MAX_ATOMS>``
  Sets the maximum number of atoms in a protein for non-Opal job submission.
  Only affects web tools.
  Default is 10000

``BUILD_PDB2PKA=False``
  Disable pkb2pka compilation.

------------------------
Web server configuration
------------------------

All the necessary files for web server installation are available with the PDB2PQR software; however, we would appreciate if users contact us before installing a publicly-accessible version of the web server so we can ensure that you are informed of PBD2PQR updates, etc.

.. note::

   These instructions are intended for systems administrators with the ability to change the behavior of their web server software and/or install software in privileged locations.

To set up a server, edit the ``build_config.py`` file and set the ``URL`` and ``PREFIX`` to appropriate values then run the usual installation procedure:

.. code-block:: bash

   $ python scons/scons.py 
   $ python scons/scons.py install

By default, the server is installed in :file:`~/pdb2pqr` and the default URL is ``http://computer_name/pdb2pqr``.

It is highly recommended that PREFIX and URL point to the same directory.
Specifying ``PREFIX=/var/www/html/pdb2pqr-test URL=http://somedomain/pdb2pqr-test`` is recommened.

If the server interface loads correctly but you cannot execute pdb2pqr by clicking the "Submit" button, make sure you have the permission to execute the :file:`pdb2pqr.cgi` file.
In particular, ensure that the access mode of :file:`pdb2pqr.cgi` allows execution by the webserver (e.g., ``chmod +x /var/www/html/pdb2pqr/pdb2pqr.cgi``).
Additionally, you may need to change the configuration of your webserver to enable CGI execution.
