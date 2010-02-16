#!/bin/bash

########################################################################
# outputfile	= actual result obtained (APBS OUTPUT)								 
# logfile		= where to place the results of the check
########################################################################

outputfile=$1
logfile=$2

docomparison(){

	if [[ $direction == "x" ]]; then
		var1=${x1}
		var2=${x2}
	fi
	
	if [[ $direction == "y" ]]; then
		var1=${y1}
		var2=${y2}
	fi
	
	if [[ $direction == "z" ]]; then
		var1=${z1}
		var2=${z2}
	fi
	
	difference=`echo $var1 $var2 | awk '{printf("%.12f",($1 - $2))}'`
	absdifference=`echo $difference | awk '{$1 = ($1 < 0) ? -$1 : $1; print}'`   # Absolute value of difference
	r=`echo $var1 $var2 $absdifference $errortol | awk '{if($1 == $2) print 1; else if($3 < $4) print 2; else print 3}'`
	
	case "$r" in
	
		1) echo "*** Comparison $comp1 for $atom1 in $direction PASSED ***"
		   echo "$outputfile: Comparison $comp1 for $atom1 in $direction PASSED ($var1)" >> $logfile
		   ;;
		2) echo "*** Comparison $comp1 for $atom1 in $direction PASSED (with rounding error - see log) ***"
		   echo "$outputfile: Comparison $comp1 for $atom1 in $direction PASSED within error ($var1; expected $var2)" >> $logfile
		   ;;
		3) echo "*** Comparison $comp1 for $atom1 in $direction FAILED ***"
	 	   echo "   APBS returned $var1"
	       echo "   Expected result is $var2 (difference of: $difference)"
		   echo "$outputfile: Comparison $comp1 for $atom1 in $direction FAILED ($var1; expected $var2)" >> $logfile
		   ;;
	esac
}

# Error tolerance as a percentage of 
errortol=0.000001

# Parse out the relevant information from the APBS output file for polar and apolar calcs
grep -A 17 'print force 1 (solv) - 2 (ref) end' $outputfile | tail -10 > polres
grep -A 17 'print APOL force 1 (asolv) end' $outputfile | tail -10 > apolres

# Merge the two files so that the calculated values and the respective know value are on the same line
# and store that in a temporary file
files="temp1 temp2"
paste polres polarforces > temp1
paste apolres apolarforces > temp2

# Now read the file line by line
for file in $files; do

	exec<$file
	while read comp1 atom1 x1 y1 z1 comp2 atom2 x2 y2 z2; do
		if [[ $x2 != "" ]]; then
			direction="x"
			docomparison $comp1 $atom1 $direction $x1 $x2
		fi
		
		if [[ $y2 != "" ]]; then
			direction="y"
			docomparison $comp1 $atom1 $direction $y1 $y2
		fi
		
		if [[ $z2 != "" ]]; then
			direction="z"
			docomparison $comp1 $atom1 $direction $z1 $z2
		fi
	done
	
done

rm -f temp* *res

exit 0

