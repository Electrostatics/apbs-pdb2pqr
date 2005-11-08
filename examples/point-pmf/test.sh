#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi


logfile=TESTRESULTS.log
nettime=0

input=( apbs ) 

results=( 1.830820736656E+01 8.906693911566E+00 5.909599752197E+00 4.430144912287E+00 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
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
  answer=`grep "Global net" ${input[i]}.out | awk '{print $5}'`
  sync

  echo ""
  echo "Energy from 1 A distance: $answer"
  if [[ ${answer} == ${results[0]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (1 A dist): PASSED ($answer)" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned $answer"
      echo "   Expected result is ${results[0]}"
      echo "           ${input[i]}.in (1 A dist): FAILED ($answer; expected ${results[0]})" >> $logfile
     
  fi
  
  # 2 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1      -1.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net" ${input[i]}.out | awk '{print $5}'`
  sync

  echo ""
  echo "Energy from 2 A distance: $answer"
  if [[ ${answer} == ${results[1]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (2 A dist): PASSED ($answer)" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned $answer"
      echo "   Expected result is ${results[1]}"
      echo "           ${input[i]}.in (2 A dist): FAILED ($answer; expected ${results[1]})" >> $logfile
  fi

  # 3 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1       0.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net" ${input[i]}.out | awk '{print $5}'`
  sync

  echo ""
  echo "Energy from 3 A distance: $answer"
  if [[ ${answer} == ${results[2]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (3 A dist): PASSED ($answer)" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned $answer"
      echo "   Expected result is ${results[2]}"
      echo "           ${input[i]}.in (3 A dist): FAILED ($answer; expected ${results[2]})" >> $logfile
  fi

  # 4 Angstrom Distance
  echo ""
  echo "ATOM      1  I   ION     1       1.000    0.000   0.000  1.00 0.00"\
     > mol2.pqr
  cat mol1.pqr mol2.pqr > complex.pqr
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=`grep "Global net" ${input[i]}.out | awk '{print $5}'`
  sync

  echo ""
  echo "Energy from 4 A distance: $answer"
  if [[ ${answer} == ${results[3]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in (4 A dist): PASSED ($answer)" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned $answer"
      echo "   Expected result is ${results[3]}"
      echo "           ${input[i]}.in (4 A dist): FAILED ($answer; expected ${results[3]})" >> $logfile
  fi


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
