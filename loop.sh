#!/bin/bash

### CALL: ./loop.sh PATH/TO/WAVE/FUNCTION/DIR/

INPUT=${1}
OUT="output.dat"

printf "time (ps)\tximp (AA)\tyimp (AA)\tzimp (AA)\tEkin (K)\tEint (K)\tESO (K)\n" > ${OUT}
for FILE in ${INPUT}density.*
do
	printf "%s\n" ${FILE}
	./get_impurity_energy.sh ${FILE} >> ${OUT}
done
