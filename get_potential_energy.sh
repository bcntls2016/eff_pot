#!/bin/bash

INPUT=${1}
ln -sf ${INPUT} djogger.dat
./4hetddft-anisotropic < input > res.dat
TIME=$(head -n 3 djogger.dat | awk '/Total evolution time/{print$5}')
EINT=$(awk '/Interaction energy/{print$5}' res.dat)
ESO=$(awk '/Spin-Orbit energy/{print$5}' res.dat)
rm djogger.dat
echo ${TIME} ${EINT} ${ESO}
