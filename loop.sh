#!/bin/bash

### CALL: ./loop.sh PATH/TO/WAVE/FUNCTION/DIR/

INPUT=${1}
OUT="output.dat"

printf "# time (ps)\tdX-HeCM (AA)\tEkin (K)\tU* (K)\t\tU+ (K)\t\tESO (K)\n" > ${OUT}
for FILE in ${INPUT}density.*
do
	printf "%s\n" ${FILE}
	./get_impurity_energy.sh ${FILE} >> ${OUT}
done
