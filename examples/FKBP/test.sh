#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi
 


logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

input=( apbs-mol apbs-smol )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: FKBP" >> $logfile
echo "Results  :" >> $logfile

# Do 1d7h-dmso first

cd 1d7h-dmso
results=( 1.500890142385E+01 1.624584726796E+01 )

# For each file in the directory, run APBS and get the value

for i in 0 1 
do
  echo "----------------------------------------"
  echo "Testing input file 1d7h-dmso/${input[i]}.in"
  echo ""

  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=( `grep "Global net" ${input[i]}.out | awk '{print $5}'` )

  echo "Global net energy: ${answer[3]}"
  sync

  fanswer=`printf "%.12f" ${answer[3]}`
  fexpected=`printf "%.12f" ${results[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`

  case "$r" in
      2)  echo "*** PASSED ***"
          echo "           1d7h-dmso/${input[i]}.in: PASSED (${answer[3]})" >> ../$logfile ;;
      1)  echo "*** PASSED (with rounding error - see log) ***"
          echo "           1d7h-dmso/${input[i]}.in: PASSED with rounding error (${answer[3]}; expected ${results[i]})" >> ../$logfile ;;
      *)  echo "*** FAILED ***"
          echo "   APBS returned ${answer[3]}"
          echo "   Expected result is ${results[i]}"
          echo "           1d7h-dmso/${input[i]}.in: FAILED (${answer[3]}; expected ${results[i]})" >> ../$logfile ;;
  esac
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

done

# Now do 1d7i-dss

cd ../1d7i-dss
results=( 1.442544510653E+01 1.545171104995E+01 )

for i in 0 1 
do
  echo "----------------------------------------"
  echo "Testing input file 1d7i-dss/${input[i]}.in"
  echo ""

  starttime=`date +%s`
  $1 ${input[i]}.in > ${input[i]}.out 
  answer=( `grep "Global net" ${input[i]}.out | awk '{print $5}'` )

  echo "Global net energy: ${answer[3]}"
  sync

  fanswer=`printf "%.12f" ${answer[3]}`
  fexpected=`printf "%.12f" ${results[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`

  case "$r" in 
      2)  echo "*** PASSED ***"
          echo "           1d7i-dss/${input[i]}.in: PASSED (${answer[3]})" >> ../$logfile ;;
      1)  echo "*** PASSED (with rounding error - see log) ***"
          echo "           1d7i-dss/${input[i]}.in: PASSED with rounding error (${answer[3]}; expected ${results[i]})" >> ../$logfile ;;
      *)  echo "*** FAILED ***"
          echo "   APBS returned ${answer[3]}"
          echo "   Expected result is ${results[i]}"
          echo "           1d7i-dss/${input[i]}.in: FAILED (${answer[3]}; expected ${results[i]})" >> ../$logfile ;;
  esac
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

done

cd ..

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
