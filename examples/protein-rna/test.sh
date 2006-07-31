#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi


logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

input=( apbs-0.050 ) 

results=(  9.941437954601E+01 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: protein-rna" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 
do
  echo "----------------------------------------"
  echo "Testing input file ${input[i]}.in"
  echo ""

 
  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out
  answer=( `grep "Global net" ${input[i]}.out | awk '{print $5}'` )
 
  sync

  fanswer=`printf "%.12f" ${answer[0]}`
  fexpected=`printf "%.12f" ${results[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`
  echo "Energy : ${answer[0]}"
  case "$r" in
      2) echo "*** PASSED ***"
         echo "           ${input[i]}.in (results): PASSED (${answer[0]})" >> $logfile ;;
      1)
         echo "*** PASSED (with rounding error - see log) ***"
         echo "           ${input[i]}.in (results): PASSED with rounding error (${answer[0]}; expected ${results[i]})" >> $logfile ;;
      *) error=`echo "scale=12;e=($fanswer - $fexpected)*100.0/$fexpected;;if(e<0)e=e*-1;e" | bc`
         ferror=`printf "%.2f" $error`        
         echo "*** FAILED ***"
         echo "   APBS returned ${answer[0]}"
         echo "   Expected result is ${results[i]} ($ferror% error)"
         echo "           ${input[i]}.in (results): FAILED (${answer[0]}; expected ${results[i]}; $ferror% error)" >> $logfile ;;
  esac
 
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
