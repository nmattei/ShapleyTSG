#!/bin/bash
set -x

#Runtime Sampler + Concorde Experiment
python3 SampleChristoCompute.py -e 1 -ml 17 -Ml 17 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/cf_exp_1_1070.csv.17only
python3 SampleChristoCompute.py -e 2 -ml 4 -Ml 17 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/cf_exp_2_1070.csv

#Error Experiments
python3 SampleChristoCompute.py -e 4 -ml 4 -Ml 10 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/cf_exp_4_1070.csv
python3 SampleChristoCompute.py -e 5 -ml 4 -Ml 10 -sa 100000 -f 1000 -sd 1 -r 1070 > ../../dat/timing_data/cf_exp_5_1070.csv


