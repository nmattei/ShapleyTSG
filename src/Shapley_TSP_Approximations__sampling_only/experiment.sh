#!/bin/bash


for i in {151..250} ; do
    rm output.$i
    touch output.$i
    
    rm output.$i.1
    touch output.$i.1
done
    
for j in {1..100} ; do 
    for i in {151..250} ; do
	
	./generate $i data.$i.$j.tsp data.$i.$j.1.tsp
	./concorde/TSP/concorde -o tour.output data.$i.$j.tsp
	./concorde/TSP/concorde -o tour.1.output data.$i.$j.1.tsp
	./mstcalculation data.$i.$j.tsp tour.output 2> /dev/null >> output.$i
	./mstcalculation data.$i.$j.1.tsp tour.1.output 2> /dev/null >> output.$i.1	
    done
done