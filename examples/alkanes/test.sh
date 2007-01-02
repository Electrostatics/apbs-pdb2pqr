#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

alkanes="2-methylbutane.pdb butane.pdb cyclohexane.pdb cyclopentane.pdb ethane.pdb hexane.pdb isobutane.pdb methane.pdb neopentane.pdb pentane.pdb propane.pdb"

results=( 1.428599279904E+01 1.220695963103E+01 1.334804666369E+01 9.557847153428E+00 9.420063118852E+00 1.622375791893E+01 1.322791252720E+01 7.793798970755E+00 1.445928389194E+01 1.445894866319E+01 1.171302556434E+01 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" >> $logfile
echo "Directory: born" >> $logfile
echo "Results  :" >> $logfile

# For each file in the directory, run APBS and get the value

i=0
for pdb in $alkanes; do

  echo "----------------------------------------"
  echo "Testing input PDB file ${pdb}"
  echo ""

  alkane=${pdb%.pdb}
  cp ${pdb} mol.pdb
  echo "Calculating for ${alkane}..."
  
  starttime=`date +%s`
  $1 apbs-apolar.in > ${alkane}.out 
  answer=`grep "Global net APOL" ${alkane}.out | awk '{print $6}'`

  echo "Global net energy: $answer"
  sync

  fanswer=`printf "%.12f" $answer`
  fexpected=`printf "%.12f" ${results[i]}`
  r=`echo "scale=12;if($fanswer>($fexpected-$vsmall) && $fanswer<($fexpected+$vsmall))r=1;if($fanswer == $fexpected)r=2;r" | bc`

  case "$r" in
      2)  echo "*** PASSED ***"
          echo "           ${alkane}.in: PASSED ($answer)" >> $logfile ;;
      1)  echo "*** PASSED (with rounding error - see log) ***"
          echo "           ${alkane}.in: PASSED with rounding error ($answer; expected ${results[i]})" >> $logfile ;;
      *)  error=`echo "scale=12;e=($fanswer - $fexpected)*100.0/$fexpected;;if(e<0)e=e*-1;e" | bc`
          ferror=`printf "%.2f" $error`
          echo "*** FAILED ***"
          echo "   APBS returned $answer"
          echo "   Expected result is ${results[i]} ($ferror% error)"
          echo "           ${alkane}.in: FAILED ($answer; expected ${results[i]}; $ferror% error)" >> $logfile ;;
  esac
  
  endtime=`date +%s`
  let elapsed=$endtime-$starttime
  let nettime=$nettime+$elapsed
  echo "Total elapsed time: $elapsed seconds"
  echo "----------------------------------------"

  i=$i+1
done

echo "Test results have been logged to ${logfile}."
echo "Total time for this directory: ${nettime} seconds."

echo "Time     : ${nettime} seconds" >> $logfile
echo "" >> $logfile
