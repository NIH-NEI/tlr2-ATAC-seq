#This scripts reads the files from genrich stdout (swarm.e), writes out the name of the sample and total peak counts

echo -e "Sample Name\tTotal Peaks Identified"

for f in *.txt
 do
	#echo $f
	samplename=$(grep "bam" $f | cut -d" " -f5 | cut -d"." -f1 | rev | cut -c2- | rev | uniq)
	numberofpeaks=$(grep "Peaks identified"  $f | cut -d" " -f3)
	echo -e "${samplename}\t${numberofpeaks}"

 done
