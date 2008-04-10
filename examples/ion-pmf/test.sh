#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

echo "This test only compares energies.  If you would like to examine forces please see README.html"

logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

coord='-3.000 -2.750 -2.500 -2.250 -2.000 -1.750 -1.50 -1.250 -1.00 -0.750 
-0.500 -0.250 0.000 0.250 0.500 0.750 1.000 1.250 1.500 1.750 2.000 
2.250 2.500 2.750 3.000'

results=( -1.494109084160E+03 -1.479876455612E+03 -1.453172513724E+03 -1.424245632015E+03 -1.395864905308E+03 -1.369209804828E+03 -1.343413508063E+03 -1.319527942713E+03 -1.297718358292E+03 -1.275175527822E+03 -1.253921668205E+03 -1.232767208837E+03 -1.211966318316E+03 -1.190891678040E+03 -1.169410509835E+03 -1.146882791843E+03 -1.125192402906E+03 -1.105192440290E+03 -1.087409226682E+03 -1.071328983874E+03 -1.056724540553E+03 -1.043946032103E+03 -1.032072744366E+03 -1.021237490030E+03 -1.011312542835E+03 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
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

  ../scripts/checkresults.sh $answer ${results[j]} $coord $logfile
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

  j=$j+1
done

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
