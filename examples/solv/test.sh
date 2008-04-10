#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi


logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

input=( apbs-mol apbs-smol ) 

methanol=( -3.624861994793E+01 -3.757593177355E+01 )
methoxide=( -3.904118544787E+02  -3.912386769429E+02 )
difference=( -3.541632345308E+02 -3.536627451694E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
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

  ../scripts/checkresults.sh ${answer[0]} ${methanol[i]} ${input[i]}.in $logfile

  # Methoxide
  
  ../scripts/checkresults.sh ${answer[1]} ${methoxide[i]} ${input[i]}.in $logfile
  
  # Difference
  
  ../scripts/checkresults.sh ${answer[2]} ${difference[i]} ${input[i]}.in $logfile
   
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
