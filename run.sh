#!/bin/sh

HOMEDIR=/home/fkirmani/Documents/CSCE798/Reproduce/svm_cpp
cd ${HOMEDIR}

make clean
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/fkirmani/Documents/CSCE798/Reproduce/svm_cpp
make
./model > features.csv

