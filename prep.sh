#!/bin/bash

INPUT=$1

cp -v ${INPUT} input.dat
sed -i "/filedenin/c\ filedenin = 'djogger.dat'" input.dat
sed -i "/mode/c\ mode = 3" input.dat
sed -i "/\&input/a\ selec_plus\t\t= '${X_PLUS}'\
 r_cutoff_plus\t= ${R_CUTOFF_PLUS}\
 umax_plus\t\t= ${UMAX_PLUS}" input.dat
