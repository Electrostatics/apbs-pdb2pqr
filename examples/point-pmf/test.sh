#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

if [[ "$3" = "ocd" ]]; then
    ocd='ocd'
else
	ocd='noocd'
fi

logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

input=( apbs ) 

results=( 1.830820736656E+01 8.906693911567E+00 5.909599752197E+00 4.430144912287E+00 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: point-pmf" >> $logfile
echo "Results  :" >> $logfile


# Initialize the fixed atom

echo "ATOM      1  I   ION     1      -3.000   0.000   0.000  1.00  0.00" \
     > mol1.pqr

# For each file in the directory, run APBS and get the value

for i in 0 
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

  starttime=`date +%s`

  # 1 Angstrom Distance

  echo "ATOM      1  I   ION     1      -2.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`
  sync

  # See if we're within a VSMALL

  ../scripts/checkresults.sh $answer ${results[0]} ${input[i]}.in $logfile $ocd
      
  # 2 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1      -1.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`
  sync

  ../scripts/checkresults.sh $answer ${results[1]} ${input[i]}.in $logfile $ocd

  # 3 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1       0.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`
  sync

  ../scripts/checkresults.sh $answer ${results[2]} ${input[i]}.in $logfile $ocd

  # 4 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1       1.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`
  sync

  ../scripts/checkresults.sh $answer ${results[3]} ${input[i]}.in $logfile $ocd

  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

done

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
