#!/bin/bash

### CALL: ./loop.sh PATH/TO/WAVE/FUNCTION/DIR/

INPUT=${1}
OUT="output.dat"

printf "time (ps)\tEkin (K)\tEint (K)\tESO (K)" >> ${OUT}
for FILE in ${INPUT}density.*
do
	./get_impurity_energy.sh ${FILE} >> ${OUT}
done