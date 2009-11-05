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

input=( apbs-mol apbs-smol )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: FKBP" >> $logfile
echo "Results  :" >> $logfile

# Do 1d7h-dmso first

cd 1d7h-dmso
results=( 1.500769108268E+01 1.624455262583E+01 )

# For each file in the directory, run APBS and get the value

for i in 0 1 
do
  echo "----------------------------------------"
  echo "Testing input file 1d7h-dmso/${input[i]}.in"
  echo ""

  starttime=`date +%s`
  ../$1 ${input[i]}.in > ${input[i]}.out 
  answer=( `grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'` )

  echo "Global net energy: ${answer[3]}"
  sync

  ../../scripts/checkresults.sh ${answer[3]} ${results[i]} 1d7h-dmso/${input[i]}.in $logfile $ocd
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

done

cat $logfile >> ../$logfile

# Now do 1d7i-dss

cd ../1d7i-dss
results=( 1.442500529301E+01 1.545150009785E+01 )

for i in 0 1 
do
  echo "----------------------------------------"
  echo "Testing input file 1d7i-dss/${input[i]}.in"
  echo ""

  starttime=`date +%s`
  ../$1 ${input[i]}.in > ${input[i]}.out 
  answer=( `grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'` )

  echo "Global net energy: ${answer[3]}"
  sync
  
  ../../scripts/checkresults.sh ${answer[3]} ${results[i]} 1d7i-dss/${input[i]}.in $logfile
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

done

cat $logfile >> ../$logfile
cd ..

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
