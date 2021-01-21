#!/bin/bash

setupATLAS

export ALRB_rootVersion=6.14.04-x86_64-slc6-gcc62-opt
lsetup root

export PYTHONPATH=$PWD/inputs/accDictionary:$PWD/scriptResonance/PythonModules:$PYTHONPATH
