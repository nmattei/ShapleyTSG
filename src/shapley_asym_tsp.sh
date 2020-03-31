#!/bin/bash

NAME=$1    # Suffix on directory name in ram filesystem.
FILENAME=$2    # Size of the game.
SAMPLES=$3 # Number of samples of game of size $SIZE
PROBLEM_SIZE=$4
PROBLEM_NAME=$5

if [ -a  /tmp/ram ] ; then 
    echo "Already have the ram file system" ; 
else 
    mkdir -p /tmp/ram
    sudo mount -t tmpfs -o size=4096M tmpfs /tmp/ram/
fi

cd /tmp/ram

CWD=`echo "/home/cgretton/Dropbox/Work2013/TSG_VRG_Research/Writing/Old_Versions/AAAI Submission/Software/Shapley_TSP_Approximations"`


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

    echo "Preparing compilation"

    mv shapley_experiment.cc shapley_experiment_OLD.cc
    echo "#define INSTALLED_PREFIX \"/tmp/ram/Concorde-$NAME/\"" > shapley_experiment.cc
    cat shapley_experiment_OLD.cc >> shapley_experiment.cc

    pwd

    echo "COMPILING"

    g++ -O6 -o one_point_matrix_experiment shapley_experiment.cc -I/home/cgretton/software/boost  -lgmpxx

    echo "DONE"

    echo "./one_point_matrix_experiment  /home/cgretton/Dropbox/Work2013/TSG_VRG_Research/Writing/Old_Versions/AAAI\ Submission/Software/$PROBLEM_NAME/Size_$PROBLEM_SIZE/$PROBLEM_NAME-$FILENAME.tsp  concorde_file_name.concorde $SAMPLES"

    #exit 0

    ./one_point_matrix_experiment  /home/cgretton/Dropbox/Work2013/TSG_VRG_Research/Writing/Old_Versions/AAAI\ Submission/Software/$PROBLEM_NAME/Size_$PROBLEM_SIZE/$PROBLEM_NAME-$FILENAME.tsp  concorde_file_name.concorde $SAMPLES

    for f in  *.data  ; do cp $f /home/cgretton/Dropbox/Work2013/TSG_VRG_Research/Writing/Old_Versions/AAAI\ Submission/Software/$PROBLEM_NAME/Size_$PROBLEM_SIZE/$PROBLEM_NAME-$FILENAME.$f ; done

    cd ..
    rm -r Concorde-$NAME
else
    echo "An error occurred while creating the directory \"Concorde-$NAME\""
fi
