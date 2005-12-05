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
  answer=( `grep "Global net" ${input[i]}.out | awk '{print $5}'` )
 
  sync

  # Methanol

  fanswer=`printf "%.12f" ${answer[0]}`
  fexpected=`printf "%.12f" ${methanol[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`
  echo "Methanol Energy : ${answer[0]}"
  case "$r" in
      2) echo "*** PASSED ***"
         echo "           ${input[i]}.in (methanol): PASSED (${answer[0]})" >> $logfile ;;
      1)
         echo "*** PASSED (with rounding error - see log) ***"
         echo "           ${input[i]}.in (methanol): PASSED with rounding error (${answer[0]}; expected ${methanol[i]})" >> $logfile ;;
      *)
         echo "*** FAILED ***"
         echo "   APBS returned ${answer[0]}"
         echo "   Expected result is ${methanol[i]}"
         echo "           ${input[i]}.in (methanol): FAILED (${answer[0]}; expected ${methanol[i]})" >> $logfile ;;
  esac

  # Methoxide

  fanswer=`printf "%.12f" ${answer[1]}`
  fexpected=`printf "%.12f" ${methoxide[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`
  echo "Methoxide Energy: ${answer[1]}"
  case "$r" in
      2) echo "*** PASSED ***"
         echo "           ${input[i]}.in (methoxide): PASSED (${answer[1]})" >> $logfile ;;
      1)
         echo "*** PASSED (with rounding error - see log) ***"
         echo "           ${input[i]}.in (methoxide): PASSED with rounding error (${answer[1]}; expected ${methoxide[i]})" >> $logfile ;;
      *)
         echo "*** FAILED ***"
         echo "   APBS returned ${answer[1]}"
         echo "   Expected result is ${methoxide[i]}"
         echo "           ${input[i]}.in (methoxide): FAILED (${answer[1]}; expected ${methoxide[i]})" >> $logfile ;;
  esac

  # Difference

  fanswer=`printf "%.12f" ${answer[2]}`
  fexpected=`printf "%.12f" ${difference[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`
  echo "Difference      : ${answer[2]}"
  case "$r" in
      2) echo "*** PASSED ***"
         echo "           ${input[i]}.in (difference): PASSED (${answer[2]})" >> $logfile ;;
      1)
         echo "*** PASSED (with rounding error - see log) ***"
         echo "           ${input[i]}.in (difference): PASSED with rounding error (${answer[2]}; expected ${difference[i]})" >> $logfile ;;
      *)
         echo "*** FAILED ***"
         echo "   APBS returned ${answer[2]}"
         echo "   Expected result is ${difference[i]}"
         echo "           ${input[i]}.in (difference): FAILED (${answer[2]}; expected ${difference[i]})" >> $logfile ;;
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
