#!/bin/bash
#echo "usage: $0 {pol:x/y/z}  {dt: in au} {emax:in eV} {npoints}"
EXE_DIR=/home/sjamnuch/Siesta/PP/
echo "usage: $0 {pol:x/y/z} {emax:in eV} {npoints}"
pd=$1
if [ -z $pd ] ;then
	pd=x;  
fi

if [[ $pd = *"x"* ]] ; then
	col=2;
elif [[ $pd = *"y"* ]] ; then
	col=3;
else 
	col=4;
fi

#dt=$2;
#if [ -z $dt ] ; then
#dt=0.08 ;
#fi

em=$2
if [ -z $em ] ; then
	em=30;
fi 

npoints=$3;
if [ -z $npoints ] ; then
	npoints=1000;
fi

em=`echo "$em/(13.605698*2)" | bc -l`;
de=`echo "$em/$npoints" | bc -l`;


dt=`grep 'deltat' dipole.vs.time | awk '{printf "%16.8f\n", $4*2}'`;

tpr=`grep 'TD.ProbeCenter' input.fdf  | awk '{print $2}'`;
if [ -z $tpr ] ; then
	tpr=$dt
else
	tpr=`echo "1.0/($tpr/(13.605698*2))" | bc -l`;
	tpr=`echo "($tpr) + ($dt)" | bc -l`;
fi

echo "TD.probecenter =" $tpr ;

#echo $col $dt $em $de
nlin=`wc dipole.vs.time | awk '{print $1}'`;
nlin=`echo "$nlin-6" | bc`;
tail -n $nlin dipole.vs.time | awk -v co="$col" '{printf "%16.8f %24.16f\n", $1/0.02418885,$co}' > j-raw.dat ;
cat j-raw.dat | awk -v delt="$tpr" '{printf "%16.8f %24.16f\n", $1-delt,$2}' > j-$pd-raw.dat ;

echo $tpr $dt
nlpre=`echo "($tpr)/($dt)" | bc -l`;
echo "nlpre =" $nlpre 
nlin=`wc j-$pd-raw.dat | awk '{print $1}'`;
echo $nlin 
nlin=`echo "scale=0;($nlin)-($nlpre)+1" | bc | awk '{printf "%ld", $1}'`;
echo $nlin

cp j-$pd-raw.dat temp.raw ;
tail -n $nlin temp.raw > j-$pd-raw.dat ;

ec=`echo "($col)+2" | bc`;
emag=`grep 'efield' dipole.vs.time | awk -v eca="$ec" '{printf "%16.8f", $eca}'`;
emag=`echo "-1*($emag)" | bc -l`;

damp=`echo "($nlin) * ($dt)/3" | bc | awk '{printf "%ld", $1}'`;

#for br in 200 150 100 75 50 ; do
for br in $damp ; do
	echo j-$pd-raw.dat $nlin $npoints $de $emag $br ;
	$EXE_DIR/epsfromj.x j-$pd-raw.dat $nlin $npoints $de $emag $br ;
	mv eps_w.dat eps_w-br$br.dat ;
done


