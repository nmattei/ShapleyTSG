#!/bin/bash
set -x
#Runtime Sampler + Concorde Experiment
python3 SampleConcordeCompute.py -e 1 -ml 16 -Ml 17 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/exp_1_1070.csv.1617
python3 SampleConcordeCompute.py -e 2 -ml 4 -Ml 17 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/exp_2_1070.csv

#Error Experiments
python3 SampleConcordeCompute.py -e 4 -ml 4 -Ml 10 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/exp_4_1070.csv
python3 SampleConcordeCompute.py -e 5 -ml 4 -Ml 10 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/exp_5_1070.csv

python3 SampleConcordeCompute.py -e 3 -ml 4 -Ml 11 -sd 1 -r 10 > ../../dat/timing_data/exp_3_11p_10.csv

