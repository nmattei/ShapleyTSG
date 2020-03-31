#!/bin/bash

X=$1

while [ $X -lt 35 ] ; do
    shapley_tsp.sh $X $X 20 12312
    
    echo $X
    X=$(($X + $2))
done
