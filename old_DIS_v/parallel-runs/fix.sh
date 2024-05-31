outfile=pwgPOWHEG+PYTHIA8-output
labels=( 'dummy1' '11' '12' '21' '22' '1H' 'H1' 'HH' )
for i in {2..8}
do
    mergedata 1 pwgPOWHEG+PYTHIA8-output-00*W$i.top -o $outfile-${labels[$i-1]}.top
done
