#!/bin/bash

INPUT=$1

cp -v ${INPUT} input.dat
sed -i "/filedenin/c\ filedenin = 'djogger.dat'" input.dat
