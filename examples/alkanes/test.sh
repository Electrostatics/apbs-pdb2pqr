#!/bin/bash

if [[ $1 = "" ]]; then
    echo "Please use \"make test\" to run the tests."
    exit
fi

logfile=TESTRESULTS.log
nettime=0
vsmall=0.000000001000

alkanes="2-methylbutane.pdb butane.pdb cyclohexane.pdb cyclopentane.pdb ethane.pdb hexane.pdb isobutane.pdb methane.pdb neopentane.pdb pentane.pdb propane.pdb"

results=( 1.439739455792E+01 1.208346456826E+01 1.354016672221E+01 9.363673200142E+00 9.422717598546E+00 1.640068943201E+01 1.323144287435E+01 7.894367190329E+00 1.449633815052E+01 1.447900211546E+01 1.192358496286E+01 )

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

forces="pentane.pdb"
for force in $forces; do

  echo "----------------------------------------"
  echo "Testing forces in ${force}"
  echo ""

  alkane=${force%.pdb}

  starttime=`date +%s`
  $1 apbs-forces.in > ${alkane}-force.out 
  
  start=`grep -n 'print APOL force 1 (solvated) end' ${alkane}-force.out | awk '{printf("%3d",$1)}'`
  stop=`grep -n 'tot   16' ${alkane}-force.out | awk '{printf("%3d",$1)}'`
  range=$(($stop-$start+1))

  grep -A ${range} 'print APOL force 1 (solvated) end' ${alkane}-force.out > answer

  difference=`diff -q answer force.result`
  
  echo 'The difference between the previous force answers and the current:'
  if [ -z "$difference" ] ; then
	echo "*** PASSED ***"
	echo "           apbs-force.in: PASSED" >> $logfile
  else
	echo "*** FAILED ***"
  	echo $difference
	echo "           apbs-force.in: FAILED $difference" >> $logfile
	echo "			 compare force.test to the standard force.result" >> $logfile
  fi
  
  cp answer force.test
  rm -f answer
  
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
