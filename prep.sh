#!/bin/bash

INPUT=$1

X_PLUS='Rb_plus_Koutselos'
R_CUTOFF_PLUS=2.d0
UMAX_PLUS=6.1105922888E+03

cp -v ${INPUT} input.dat
sed -i "/filedenin/c\ filedenin\t\t= 'djogger.dat'" input.dat
sed -i "/mode/c\ mode\t\t\t= 3" input.dat
sed -i "/\&input/a\ selec_plus\t\t= '${X_PLUS}'\n \
r_cutoff_plus\t\t= ${R_CUTOFF_PLUS}\n \
umax_plus\t\t= ${UMAX_PLUS}" input.dat
