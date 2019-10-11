# Initial UQ tools for APBS
The initial tools to perform UQ with APBS have been developed using Matlab. This Matlab code, along with a Python wrapper can be downloaded using the instructions from the following section.

## Download Instructions
You will need to have Python 3.5 or higher interpreter installed in your system. Click [here](https://www.python.org/) to download the one for your system.

To get Matlab and related products follow [this link](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml). Once Matlab is installed on your system you will need to install Matlab's engine API for python. You may need admin privileges. To install the engine follow the steps below. More detailed instructions can be found [here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).

* From a command prompt window:

  `$ cd <matlab_root_dir>/extern/engines/python`

  `$ python setup.py install`

## Running the UQ tools for APBS

The tool can be run from a command prompt or any other environment that will allow you to run python and pass command line arguments. To run enter at the command prompt:

`python dsrar.py [-h] --fileName FILENAME --pqrName PQRNAME
                [--polyOrder POLYORDER] [--Nrandomdim NRANDOMDIM]
                [--Npart NPART] [--Nperpart NPERPART] [--startN STARTN]
                [--stopN STOPN] [--stepS STEPS] [--procs PROCS]`

For a description of the flags and options enter:

`python dsrar.py --help`
