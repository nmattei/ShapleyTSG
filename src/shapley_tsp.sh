#!/bin/bash

NAME=$1    # Suffix on directory name in ram filesystem.
SIZE=$2    # Size of the game (players).
SAMPLES=$3 # Number of samples of game of size $SIZE
SEED=$4

if [ -a  /tmp/ram ] ; then 
    echo "Already have the ram file system" ; 
else 
    mkdir -p /tmp/ram
    sudo mount -t tmpfs -o size=2048M tmpfs /tmp/ram/
fi

cd /tmp/ram

CWD=`echo "/home/cgretton/Dropbox/Work2013/TSG_VRG_Research/Writing/Old_Versions/AAAI Submission/Software/Shapley_TSP_Approximations/"`
#CWD="~/Dropbox/Work2013/Software/Shapley_TSP_Approximations/"

if [ -a Concorde-$NAME  ] ; then 
    echo "Directory \"Concorde-$NAME\" already exists, exiting!"
else 
    echo "Executing: cp -r "$CWD"  Concorde-$NAME"
    cp -r "$CWD"  Concorde-$NAME
    #DIR_NAME=`echo $CWD | sed s/"\/"/" "/g | awk '{print $NF}'`
    #mv $DIR_NAME Concorde-$NAME
fi

if [ -a Concorde-$NAME  ] ; then 
    cd Concorde-$NAME

    mv shapley_experiment.cc shapley_experiment_OLD.cc
    echo "#define INSTALLED_PREFIX \"/tmp/ram/Concorde-$NAME/\"" > shapley_experiment.cc
    cat shapley_experiment_OLD.cc >> shapley_experiment.cc

    echo "g++ -O6 -o shapley_experiment shapley_experiment.cc -I/home/cgretton/software/boost  -lgmpxx  -lgmp"
    g++ -O6 -o shapley_experiment shapley_experiment.cc -I/home/cgretton/software/boost  -lgmpxx -lgmp

    # 7 JAN 2014 #
    #./shapley_experiment $SIZE $SIZE $SAMPLES concorde_file_name.concorde 1000 1000 2000 0 12312
    # 9 JAN 2014 #./shapley_experiment $SIZE $SIZE $SAMPLES concorde_file_name.concorde 1000 1000 2000 0 21321
    ./shapley_experiment $SIZE $SIZE $SAMPLES concorde_file_name.concorde 100 300000 300000 0 $SEED

    for f in *.data ; do cp $f /home/cgretton/Dropbox/Data/$f-$NAME ; done
    for f in *.tsp ; do cp $f /home/cgretton/Collection/ ; done
    
    cd ..
    ##rm -r /tmp/ram/Concorde-$NAME
else
    echo "An error occurred while creating the directory \"Concorde-$NAME\""
fi
