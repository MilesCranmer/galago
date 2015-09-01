#!/bin/bash
for x in $(ls *.bin) 
	do
	y=$(echo $x | cut -f3 -d"_" | cut -c 1-3)
	if [ $y -lt 270 ]
		then 
			mv $x weirdos/
	fi
done

