#!/bin/bash

SIZE=$1 ## i.e. 36
SAMPLES=$2


if [ $SAMPLES -eq 1 ] ; then
    for i in {1..16} ; do
	mkdir at50_S1.$i
	cd at50_S1.$i
	../shapley_experiment $SIZE $SIZE 1 concorde_file_name.concorde 100 1000 0 $i  &
	cd ..
    done

elif [ $SAMPLES -eq 2 ] ; then
    for i in {17..32} ; do
	mkdir at50_S2.$i
	cd at50_S2.$i
	../shapley_experiment $SIZE $SIZE 1 concorde_file_name.concorde 100 1000 0 $i  &
	cd ..
    done
else
    for i in {16..20} ; do
	mkdir at50_S3.$i
	cd at50_S3.$i
	../shapley_experiment $SIZE $SIZE 1 concorde_file_name.concorde 100 1000 0 $i &
	cd ..
    done
fi

