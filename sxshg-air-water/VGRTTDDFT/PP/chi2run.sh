#!/bin/bash


for i in [0-9]*
do
	cd $i
for j in *E*/
do
	cd $j
	awk '{print $1,sqrt($2**2+$3**2)}' ftout* > jw.dat
	w=`sort -g -k2 jw.dat | tail -n 1`
	e=`echo $i | awk '{print $1*2}'`
	tw=`sort -g -k2 jw.dat | awk -v e=$e '{if ( $1 > e- 20  && $1 < e+20) print $0}' | tail -n 1`
	echo $j $w $tw
	cd ../
done > CHI2_raw.dat
sort -g -k1 CHI2_raw.dat > CHI2.dat
sed -i "s/\///g" CHI2.dat
rm CHI2_raw.dat
cd ../
done
