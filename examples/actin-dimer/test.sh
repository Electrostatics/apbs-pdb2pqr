#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

if [[ $2 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

logfile=TESTRESULTS.log
nettime=0

input=( apbs-mol-auto apbs-smol-auto apbs-spl2-auto apbs-mol-parallel apbs-smol-parallel apbs-spl2-parallel )

results=( 1.131989619035E+02 1.129865798917E+02 1.445817673653E+01 9.863148054323E+01 1.173921272889E+02 2.060009725426E+01 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: actin-dimer" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 1 2 3 4 5
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

  procs=`grep -m 1 pdime ${input[i]}.in | awk '{print $2*$3*$4}'`
  if [[  $procs != "" ]]; then
      # Split the files

      echo "Splitting the input file into $procs separate files using the"
      echo "inputgen.py utility..."
      echo ""

      python $2 --split ${input[i]}.in
      
      j=0
      answer=0
      starttime=`date +%s`
      while [ $j -lt $procs ]
      do
        $1 ${input[i]}-PE$j.in > ${input[i]}-PE$j.out
        a=`grep "Global net" ${input[i]}-PE$j.out | awk '{print $5}'`
        echo "Processor $j result: $a"
        conv=`printf "%.12f" $a`
        answer=`echo "scale=12; $answer+$conv" | bc`
        echo ""
        let j=$j+1
      done
      answer=`printf "%.12E" $answer`

  else
      starttime=`date +%s`
      $1 ${input[i]}.in > ${input[i]}.out 
      answer=`grep "Global net" ${input[i]}.out | awk '{print $5}'`
  fi


  echo "Global net energy: $answer"
  sync
  if [[ $answer = ${results[i]} ]]; then
      echo "*** PASSED ***"
      echo "           ${input[i]}.in: PASSED ($answer)" >> $logfile
  else
      echo "*** FAILED ***"
      echo "   APBS returned $answer"
      echo "   Expected result is ${results[i]}"
      echo "           ${input[i]}.in: FAILED ($answer; expected ${results[i]})" >> $logfile
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
