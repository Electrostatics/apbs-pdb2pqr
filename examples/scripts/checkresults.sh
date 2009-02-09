#!/bin/bash

########################################################################
# result	= actual result obtained								 
# expected	= results based on default test config
# testfile	= name of the .in file used in this run
########################################################################

result=$1
expected=$2
testfile=$3
logfile=$4

# Error tolerance as a percentage of 
errortol=0.0001

error=`echo $result $expected | awk '{printf("%.12f",(($1 - $2)*100)/$2)}'`
relativeError=`echo $error | awk '{$1 = ($1 < 0) ? -$1 : $1; print}'`

r=`echo $result $expected $relativeError $errortol | awk '{if($1 == $2) print 1; else if($3 < $4) print 2; else print 3}'`

case "$r" in

	1) echo "*** PASSED ***"
	   echo "$testfile: PASSED ($result)" >> $logfile
	   ;;
	2) echo "*** PASSED (with rounding error - see log) ***"
	   echo "$testfile: PASSED within error ($result; expected $expected; $relativeError% error)" >> $logfile
	   ;;
	3) echo "*** FAILED ***"
 	   echo "   APBS returned $result"
           echo "   Expected result is $expected ($error% error)"
	   echo "$testfile: FAILED ($result; expected $expected; $relativeError% error)" >> $logfile
	   ;;
esac

exit 0

