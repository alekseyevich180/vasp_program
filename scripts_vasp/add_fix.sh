#!/bin/bash
c=$1
cp POSCAR POSCAR_orig
rm POSCAR

selective=0

str1="Selective"

while read line
do
str2=`eval echo $line | awk '{print $1}'`
if [ "$str1" == "$str2" ]
then
    selective=1
else
    :
fi
    done < POSCAR_orig

if [ $selective -eq 1 ]; then
awk 'NR <= 9 {print $0}' POSCAR_orig >> POSCAR
awk 'NR >= 10 {if($3 < '$c'){print "     " $1 "     " $2 "     " $3 "     F F F     "}else{print "     " $1 "     " $2 "     " $3 "     T T T     "}}' POSCAR_orig >> POSCAR
else
awk 'NR <= 7 {print $0}' POSCAR_orig >> POSCAR
echo 'Selective dynamics' >> POSCAR
awk 'NR == 8 {print $0}' POSCAR_orig >> POSCAR
awk 'NR >= 9 {if($3 < '$c'){print "     " $1 "     " $2 "     " $3 "     F F F     "}else{print "     " $1 "     " $2 "     " $3 "     T T T     "}}' POSCAR_orig >> POSCAR
fi

sed -i -e "s// /g" POSCAR
