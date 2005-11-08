#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi


logfile=TESTRESULTS.log
nettime=0

input=( apbs-mol apbs-smol apbs-spl2 ) 

methanol=( -3.624861994793E+01 -3.757593177355E+01  -6.434717255302E+01 )
methoxide=( -3.904118544787E+02  -3.912386769429E+02 -4.153175490634E+02 )
difference=( -3.541632345308E+02 -3.536627451694E+02 -3.509703765104E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: solv" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 1 2 
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

 
  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out
  answer=( `grep "Global net" ${input[i]}.out | awk '{print $5}'` )
 
  sync

  echo "Methanol Energy : ${answer[0]}"
  if [[ ${answer[0]} == ${methanol[i]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (methanol): PASSED" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned ${answer[0]}"
      echo "   Expected result is ${methanol[i]}"
      echo "           ${input[i]}.in (methanol): FAILED (${answer[0]}; expected ${methanol[i]})" >> $logfile
  fi

  echo "Methoxide Energy: ${answer[1]}"
  if [[ ${answer[1]} == ${methoxide[i]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (methoxide): PASSED" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned ${answer[1]}"
      echo "   Expected result is ${methoxide[i]}"
      echo "           ${input[i]}.in (methoxide): FAILED (${answer[1]}; expected ${methoxide[i]})" >> $logfile
  fi

  echo "Difference      : ${answer[2]}"
  if [[ ${answer[2]} == ${difference[i]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (difference): PASSED" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned ${answer[2]}"
      echo "   Expected result is ${difference[i]}"
      echo "           ${input[i]}.in (difference): FAILED (${answer[2]}; expected ${difference[i]})" >> $logfile
  fi
  
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
