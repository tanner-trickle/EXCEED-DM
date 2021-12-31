#!/bin/bash

EXDM_PATH="./build/exdm"

for dir in ./examples/*/ ; do


    if [ $dir != './examples/data/' ]
    then

        mpirun -np $1 $EXDM_PATH "${dir}input.txt"

    fi

done

