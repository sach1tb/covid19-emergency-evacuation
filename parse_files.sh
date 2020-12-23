# usage parse_files.sh <parameter file> <number of runs>
# eg. parse_files.sh sd_u1s1.par 10
set -v
for ii in `seq $2`
do
	tmpfile=$(printf "%s_%.2d.par" $(basename $1 .par) $ii)
	echo $tmpfile
	# id=$(echo "$tmpfile" | sed "s/\.par//g")
	id=$(basename "$tmpfile" .par)
	echo $id
	sed "s/RndSeed.*123/RndSeed   $ii/g" $1 > $tmpfile
	./sd 0 $tmpfile
        mv $tmpfile dump/
	mv sd.dat dump/$id.dat
	mv sd.dat2 dump/$id.dat2
	mv sd.dat3 dump/$id.dat3
done
