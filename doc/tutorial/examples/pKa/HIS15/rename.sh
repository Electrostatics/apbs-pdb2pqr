#!/bin/bash

for i in *.in; do
	new=`echo $i | sed -e "s/GLH35/HSP15/g" | sed -e "s/GLU35/HIS15/g"`
	cat $i | sed -e "s/GLH35/HSP15/g" | sed -e "s/GLU35/HIS15/g" > ${new}
	rm $i
	echo $i $new
done
