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

echo "This test only compares energies.  If you would like to examine forces please see README.html"

logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

coord='1.000'

results=( -1.125192402906E+03 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: ion-pmf" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

j=0
for i in $coord
do
  echo "----------------------------------------"
  echo "Testing ion at $i"
  echo ""

  # TO TRANSLATE IN THE X DIRECTION:
  echo "ATOM      1  ION   ION     1      -3.000   0.000   0.000  1.00 2.00" > complex.pdb
  echo "ATOM      1  ION   ION     1      $i    0.000   0.000  1.00 2.00" >> complex.pdb
  
  starttime=`date +%s`
  $1 apbs.in > OUTPUT_${i}.out 
  answer=`grep "Global net ELEC" OUTPUT_${i}.out | awk '{print $6}'`

  echo "Global net energy: $answer"
  sync

  ../scripts/checkresults.sh $answer ${results[j]} apbs${i} $logfile $ocd
  ../scripts/checkforces.sh OUTPUT_${i}.out $logfile
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

  let j=$j+1
done

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile


