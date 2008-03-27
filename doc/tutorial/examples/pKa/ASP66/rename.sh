#!/bin/bash

for i in *.in; do
	new=`echo $i | sed -e "s/GLH35/ASH66/g" | sed -e "s/GLU35/ASP66/g"`
	cat $i | sed -e "s/GLH35/ASH66/g" | sed -e "s/GLU35/ASP66/g" > ${new}
	rm $i
	echo $i $new
done
