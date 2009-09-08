#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

if [[ $2 = "" ]]; then
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

input=( apbs-smpbe-24dup apbs-smpbe-ion )

results=( 8.404388689915E+03 8.141539945656E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: smpbe" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

for i in 0 1
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
        a=`grep "Global net ELEC" ${input[i]}-PE$j.out | awk '{print $6}'`
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
      answer=`grep "Global net ELEC" ${input[i]}.out | awk '{print $6}'`
  fi
 
  echo "Global net energy: $answer"
  sync
  
  ../scripts/checkresults.sh $answer ${results[i]} ${input[i]}.in $logfile $ocd
    
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
