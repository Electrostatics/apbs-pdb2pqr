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

acetic=( -2.267881997628E+01 -2.233050451129E+01 ) 
acetate=( -1.997462580204E+02 -1.984883191396E+02 ) 
proton=( -2.974598331751E+02 -2.959668653531E+02 )
ionization=( -4.745272868358E+02 -4.721247084138E+02 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: ionize" >> $logfile
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

  # Acetic Acid
  
  echo "Acetic Acid Energy: ${answer[0]}"
  ../scripts/checkresults.sh ${answer[0]} ${acetic[i]} ${input[i]}.in $logfile $ocd

  # Acetate

  echo "Acetate Energy    : ${answer[1]}"
  ../scripts/checkresults.sh ${answer[1]} ${acetate[i]} ${input[i]}.in $logfile $ocd
  
  # Proton

  echo "Proton Energy     : ${answer[2]}"
  ../scripts/checkresults.sh ${answer[2]} ${proton[i]} ${input[i]}.in $logfile $ocd

  # Ionization

  echo "Ionization Energy : ${answer[3]}"
  ../scripts/checkresults.sh ${answer[3]} ${ionization[i]} ${input[i]}.in $logfile $ocd

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
