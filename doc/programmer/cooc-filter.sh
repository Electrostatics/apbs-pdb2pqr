#!/bin/bash

# Check invocation
if test $# != 1; then 
    echo " Got $# arguments, need 1."
    echo "\n Usage:  c2c++ <file>"
fi

# Get the file name
filename=$1

cat $filename | sed "s/VEXTERNC //g"
