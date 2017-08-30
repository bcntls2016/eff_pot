#!/bin/bash

### CALL: ./get_impurity_energy.sh path/to/wave-function

INPUT=${1}
ln -sf ${INPUT} djogger.dat
./4hetddft-anisotropic < input > res.dat
TIME=$(head -n 3 djogger.dat | awk '/Total evolution time/{print$5}')
EKIN=$(awk '/Kinetic energy \(X/{print$5}' res.dat)
EINT=$(awk '/Interaction energy/{print$5}' res.dat)
ESO=$(awk '/Spin-Orbit energy/{print$5}' res.dat)
rm djogger.dat
printf "%.1f\t%f\t%f\t%f" ${TIME} ${EKIN} ${EINT} ${ESO}
