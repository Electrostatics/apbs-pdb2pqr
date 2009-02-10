#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi



logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

input=( apbs-mol-vdw apbs-smol-vdw apbs-mol-surf apbs-smol-surf )

results=( 8.085852882999E+00 2.096282333147E+01 1.192607452865E+02 1.088773806893E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: pka-lig" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 1 2 3 
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`

  echo "Global net energy: $answer"
  sync

  ../scripts/checkresults.sh $answer ${results[i]} ${input[i]}.in $logfile
  
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
