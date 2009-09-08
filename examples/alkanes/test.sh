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

alkanes="2-methylbutane.pdb butane.pdb cyclohexane.pdb cyclopentane.pdb ethane.pdb hexane.pdb isobutane.pdb methane.pdb neopentane.pdb pentane.pdb propane.pdb"

results=( 1.439739455792E+01 1.208346456826E+01 1.354016672221E+01 9.363673200142E+00 9.422717598546E+00 1.640068943201E+01 1.323144287435E+01 7.894367190329E+00 1.449633815052E+01 1.447900211546E+01 1.192358496286E+01 )

# Initialize the results file

date=`date`
echo "Date     : ${date}" > $logfile
echo "Directory: alkanes" >> $logfile
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

  ../scripts/checkresults.sh $answer ${results[i]} ${alkane}.in $logfile $ocd
  
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
