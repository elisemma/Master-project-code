#!/bin/bash

FILE=46Ti_coh.dat
BIN=coh
DIR=output/

OUTBEG=out_coh_
OUTEXT=_MeV.dat

# get B2 values from:
# https://t2.lanl.gov/nis/molleretal/publications/ADNDT-FRDM2012.pdf

# Remove old output CoH* files
rm -f CoH*
# Remove old calculations from output folder
rm -r $DIR/out_coh_*

# adjust sequence accordingly
for E in $(seq 1 2.5 55)
do
    echo "Calculating $E -MeV incident deuterons";
#    touch $(printf %03.0f $E)$OUTEXT;
#    $BIN -f -p 19 -e $E < $FILE > $DIR$OUTBEG$E$OUTEXT;
    $BIN -f -p 19 -e $E < $FILE > $DIR$OUTBEG$(printf %03.0f $E)$OUTEXT;
done
