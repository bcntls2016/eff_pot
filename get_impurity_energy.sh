#!/bin/bash

### CALL: ./get_impurity_energy.sh PATH/TO/WAVE/FUNCTION/FILE

INPUT=${1}

ln -sf ${INPUT} djogger.dat
./4hetddft-anisotropic < input.dat > res.dat
TIME=$(head -n 3 djogger.dat | awk '/Total evolution time/{print$5}')
XIMP=$(head djogger.dat | awk 'NR==9{print$1}')
YIMP=$(head djogger.dat | awk 'NR==9{print$2}')
ZIMP=$(head djogger.dat | awk 'NR==9{print$3}')
EKIN=$(awk '/Kinetic energy \(X/{print$5}' res.dat)
EINT=$(awk '/Interaction energy \(X\*/{print$5}' res.dat)
EINTPLUS=$(awk '/Interaction energy \(X\+/{print$5}' res.dat)
ESO=$(awk '/Spin-Orbit energy/{print$5}' res.dat)
rm djogger.dat res.dat DFT4He3d.namelist.read
printf "%.1f\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" ${TIME} ${XIMP} ${YIMP} ${ZIMP} ${EKIN} ${EINT} ${EINTPLUS} ${ESO}
