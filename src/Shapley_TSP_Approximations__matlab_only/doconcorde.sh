#!/bin/bash

echo "./concorde/TSP/concorde $1"

X=`./concorde/TSP/concorde $1 2> /dev/null | grep "Optimal Solution:" | awk '{print $3}'`



exit $X
