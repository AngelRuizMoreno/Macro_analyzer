#!/bin/bash
for f in *.sdf;
do
	echo $f 
	t=$(date +"%T")
	echo "start: $t"
	./Mcnf -w0 -c3 -v2 $f
	l=$(date +"%T")
	echo "finish: $l"
done
