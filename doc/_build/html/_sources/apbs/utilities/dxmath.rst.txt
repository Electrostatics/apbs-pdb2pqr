.. _dxmath:

dxmath
======

Performs simple arithmetic operations with Cartesian grid data.  
This program takes as input a file with operations specified in a stack-based (RPN) manner.
For example, a command file which adds grid1 and grid2, multiplies the result by 5.3, adds grid4, subtracts 99.3 from the whole thing, and writes the result on grid5 would have the form:

.. code-block:: mathematica
   
   grid1
   grid2 +
   5.3 *
   grid4 +
   99.3 -
   grid5 =

The file names, scalar values, and operations must be separated by tabs, line breaks, or white space.
Comments can be included between the character # and a new line (in the usual shell script fashion).
Found in :file:`tools/mesh`
