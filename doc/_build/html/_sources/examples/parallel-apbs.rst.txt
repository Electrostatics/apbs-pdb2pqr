Parallel APBS execution for large calculations
==============================================

=============
Why parallel?
=============

APBS finite difference multigrid calculations require approximately 200 B memory per grid point.
These memory requirements can be distributed in two ways during a calculation:

* APBS calculations can be performed in parallel across multiple processors (hopefully, sharing distributed memory!). This functionality is provided by using the :ref:`mgpara` keyword.

* APBS calculations can be broken into a series of smaller, asynchronous runs which (individually) require less memory. This functionality is provided by using both the :ref:`mgpara` and :ref:`async` keywords. 

=================================
Synchronous parallel calculations
=================================

The actin dimer example provided with the APBS distribution :file:`examples/actin-dimer/` is a fairly large system that can often require too much memory for some systems. 
This example will use the actin dimer complex PQR file (:file:`complex.pqr`) to illustrate parallel focusing.

We're going to use an 8-processor parallel calculation to write out the electrostatic potential map for this complex.
Each processor will solve a portion of the overall problem using the parallel focusing method on a 973 mesh with 20% overlap between meshes for neighboring processors.
An example input file for this calculation might look like:

.. code-block:: bash
   
   read
     mol pqr complex.pqr
   end
   elec name complex
     mg-para
     ofrac 0.1
     pdime 2 2 2
     dime 97 97 97
     fglen 150 115 160
     cglen 156 121 162
     cgcent mol 1
     fgcent mol 1
     mol 1
     npbe
     bcfl sdh
     ion 1 0.150 2.0
     ion -1 0.150 2.0 
     pdie 2.0
     sdie 78.54
     srfm mol
     chgm spl0
     srad 1.4
     swin 0.3
     sdens 10.0
     temp 298.15
     calcenergy total
     calcforce no
     write pot dx pot
   end
   quit

where the ":ref:`pdime` 2 2 2" statement specifies the 8-processor array dimensions, the ":ref:`ofrac` 0.1" statement specifies the 20% overlap between processor calculations, and the ":ref:`dime` 97 97 97` statement specifies the size of each processor's calculation.
The ":ref:`write` pot dx potential" instructs APBS to write out OpenDX-format maps of the potential to 8 files :file:`potential-{#}.dx`, where *#* is the number of the particular processor.

An MPI-compiled version of APBS can be used with this input file to run 8 parallel focusing calculations, with each calculation generating fine-scale solutions on a different region of the (:ref:`fglen`) problem domain.
Note that 8 separate OpenDX files are written by the 8 processors used to perform the calculation.
Writing separate OpenDX< files allows us to avoid communication in the parallel run and keeps individual file sizes (relatively) small.
Additionally, if a user is interested in a specific portion of the problem domain, only a few files are needed to get local potential information.
However, most users are interested in global potentials. 
APBS provides the :ref:`mergedx` program to reassemble the separate OpenDX files into a single file.
`mergedx` is a simple program that allows users to combine several OpenDX files from a parallel focusing calculation into a single map.
This map can be down-sampled from the original resolution to provide coarser datasets for fast visualization, etc.
For example, the command

.. code-block:: bash
   
   $ mergedx 65 65 65 pot0.dx pot1.dx pot2.dx pot3.dx pot4.dx pot5.dx pot6.dx pot7.dx

will generate a file :file:`gridmerged.dx` which has downsampled the much larger dataset contained in the 8 OpenDX files into a 65\ :sup:`3` file which would be suitable for rough visualization.
An example of mergedx output visualization is shown in the attached figure.
Note that downsampling isn't necessary -- and often isn't desirable for high quality visualization or quantitative analysis.

.. image:: /media/actin_dimer-iso_trans.jpg

==================================
Asynchronous parallel calculations
==================================

The steps described in the previous section can also be performed for systems or binaries which are not equipped for parallel calculations via MPI.
In particular, you can add the statement ":ref:`async` *n*" to the ELEC :ref:`mgpara` section of the APBS input file to make the single-processor calculation masquerade as processor *n* of a parallel calculation.

Scalar maps from asynchronous APBS calculations can be combined using the mergedx program as described above.
Currently, energies and forces from asynchronous APBS calculations need to merged manually (e.g., summed) from the individual asynchronous calculation output.
This can be accomplished by simple shell scripts.

