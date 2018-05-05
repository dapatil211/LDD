#!/bin/bash

GRAPH=$1
OUT=$2

BETAS=(".002" ".005" ".01" ".02" ".05" ".1" ".2")
#BETAS=(".002" ".005")

#NUMTHREADS=("1" "2" "4" "8" "16" "36")
NUMTHREADS=("1" "2" "4")

for BETA in ${BETAS[*]}
do
./LDDSerial -e $BETA -s $GRAPH > $OUT"_serial_"$BETA
done

for BETA in ${BETAS[*]}
do
for NUMTHREAD in ${NUMTHREADS[*]}
do
./LDD -p $NUMTHREAD -e $BETA -s $GRAPH > $OUT"_par_"$NUMTHREAD"_"$BETA
done
done
echo DONE >&2