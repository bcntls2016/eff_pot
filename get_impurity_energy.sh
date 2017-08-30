#!/bin/bash

INPUT=${1}
ln -sf ${INPUT} djogger.dat
./4hetddft-anisotropic < input > res.dat
TIME=$(head -n 3 djogger.dat | awk '/Total evolution time/{print$5}')
EKIN=$(awk '/Kinetic energy \(X/{print$5}' res.dat)
EINT=$(awk '/Interaction energy/{print$5}' res.dat)
ESO=$(awk '/Spin-Orbit energy/{print$5}' res.dat)
rm djogger.dat
echo ${TIME} ${EKIN} ${EINT} ${ESO}
printf "%d\t%f\t%d" ${TIME} ${EKIN} ${EINT} ${ESO}
