#!/bin/bash

./4hetddft-anisotropic < input > res.dat
TIME=$(head -n 3 djogger.dat | awk '/Total evolution time/{print$5}')
EINT=$(awk '/Interaction energy/{print$5}' res.dat)
echo ${TIME} ${EINT}
