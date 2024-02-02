export input=$1  # EXPECT ID\tLISTOFREGIONS
export genes=$2  # EXPECT CHROM\tSTART\tEND\tID
export output=$3
export ScriptDir=$4

perl $ScriptDir/vlookup.pl $input 0 $genes 3 0 $input.wChroms
if [ -s $output ]
then
    rm $output
fi
touch $output
export counter=1
while read prom
do
  export id=$(echo $prom | awk -F' ' '{print $1}')
  export chrom=$(echo $prom | awk -F' ' '{print $3"_"}')
  export editedRegions=$(echo $prom | awk -F' ' '{print $2}' | sed s/\,/\\n/g | awk -F\\t -v chrom=$chrom '{if($1 ~ chrom) print $0}' | sed ':a;N;$!ba;s/\n/,/g')

  echo -e "$id\t$editedRegions" >> $output
  #### COUNT PROCESSED GENES SO FAR
    counter=$((counter+1))
    if ! (($counter % 1000)); then
        echo "$counter promoters processed"
    fi
done < $input.wChroms
