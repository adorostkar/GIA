#!/bin/sh
# get files with gnuplot ending
FILES=*.gnuplot

# create the title array
i=1;
for f in $FILES
do
	NA="${f/surface/nue=}"
	NA="${NA/_/}"
	NA="${NA/_/.}"
	NA="${NA/.gnuplot/}"
	NAMES[$i]=$NA
	let i=i+1
done

#pass everything to gnuplot
gnuplot << EOF
reset
set title 'Vertical displacement'
set xlabel 'X values(Km)'
set ylabel 'Vertical displacemnets(m)'
set key right bottom

filenames=system("ls -1 *gnuplot")
titl="${NAMES[@]}"
count=words(filenames)
files(i)=word(filenames,i)
titles(i)=word(titl,i)

plot for [i=1:count] files(i) using (\$1*10000):(\$3) with lines title titles(i)

EOF

