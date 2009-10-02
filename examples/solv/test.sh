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

methanol=( -3.624863445503E+01 -3.757593797629E+01 )
methoxide=( -3.904121297757E+02  -3.912388198513E+02 )
difference=( -3.541634953207E+02 -3.536628818750E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: solv" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 1 
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

 
  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out
  answer=( `grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'` )
 
  sync

  # Methanol

  ../scripts/checkresults.sh ${answer[0]} ${methanol[i]} ${input[i]}.in $logfile $ocd

  # Methoxide
  
  ../scripts/checkresults.sh ${answer[1]} ${methoxide[i]} ${input[i]}.in $logfile $ocd
  
  # Difference
  
  ../scripts/checkresults.sh ${answer[2]} ${difference[i]} ${input[i]}.in $logfile $ocd
   
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

  # Put the results of all the tests in TEST-RESULTS.log

done

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
