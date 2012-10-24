#!/bin/bash

rm -rf ../bin/* ../tools/bin/* CMakeFiles CMakeCache.txt cmake_install.cmake Makefile src install_manifest.txt apbs.cbp maloc-prefix tools fontconfig ../tests/*.pyc ../tests/io.mc ../tests/*~ ../tests/test.log doc ../doc/programmer/html

find .. -name *~ | xargs rm

