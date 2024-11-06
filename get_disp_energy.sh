#!/bin/bash
B=`awk 'END {print $3}' fe.dat`
C=`cat OUTCAR | grep Edisp | tail -1 | awk '{print $3}'`
echo $i $B  $C >> energy.txt
