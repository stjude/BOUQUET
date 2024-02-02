
module load bedtools/2.30.0

echo -e "\n\nFILES TO MERGE:\n$@"

# multijoin - join multiple files

join_rec() {
    if [ $# -eq 1 ]; then
        join -t$'\t' -a1 -a2 - <(sort -k1,1 -t$'\t' "$1")
    else
        f=$1; shift
        join -t$'\t' -a1 -a2 - <(sort -k1,1 -t$'\t' "$f") | join_rec "$@"
    fi
}


WD=$(dirname $1)
if [ $# -eq 1 ]; then
    sort -t$'\t' -k1,1 "$1"|sed -e 's/\t/,/g'|sed -e 's/,/\t/1' > $WD/temp
elif [ $# -eq 2 ]; then
    join -t$'\t' -a1 -a2 <(sort -t$'\t' -k1,1 "$1") <(sort -t$'\t' -k1,1 "$2")|sed -e 's/\t/,/g'|sed -e 's/,/\t/1' > $WD/temp
else
    f1=$1; f2=$2; shift 2
    join -t$'\t' -a1 -a2 <(sort -t$'\t' -k1,1 "$f1") <(sort -t$'\t' -k1,1 "$f2") | join_rec "$@"|sed -e 's/\t/,/g'|sed -e 's/,/\t/1' > $WD/temp
fi


: << 'END'

myArray=( "$@" )
export iterator=1

sort -k1,1n ${myArray[0]} > temp

for arg in "${myArray[@]}"
do
  export files=$(echo -e "$files\t$arg")
  sort -k1,1n $arg > temp$iterator

  paste temp temp$iterator | awk -F\\t '{print $1 "\t" $2 "," $4}' > tempX
  mv tempX temp

  iterator=$((iterator+1))
done

exit
END

cd $WD

if [ -s RegionsCollapsed.txt ]
then
  rm RegionsCollapsed.txt
fi
touch RegionsCollapsed.txt

export counter=1
while read promLine
do
  # tab was tranformed to space by echo function
  export promoter=$(echo $promLine | awk -F' ' '{print $1}')
  echo -e "$promoter" > thisProm.txt
#   echo $promLine
  echo $promLine | awk -F' ' '{print $2}' > theseRegions.txt
  #remove empty lines by sed, add by Jie, 11/02/2020
  sed s/\,/\\n/g theseRegions.txt |sed '/^[[:space:]]*$/d'| sed s/\_/\\t/g | sort -k1,1 -k2,2n > theseRegions.bed
  #sed s/\,/\\n/g theseRegions.txt | sed s/\_/\\t/g | awk -F\\t '{if($1 != "") print $1 "\t" $2 "\t" $3}' | sort -k1,1 -k2,2n > theseRegions.bed
  #do not merge neighbour bins togenter, keep location to be individual bin (5K)
  #mergeBed -i theseRegions.bed > theseRegionsCollapsed.bed
  #sed s/\\t/\_/g theseRegionsCollapsed.bed | sed ':a;N;$!ba;s/\n/,/g' > theseRegionsCollapsed.txt
  sed s/\\t/\_/g theseRegions.bed | sed ':a;N;$!ba;s/\n/,/g' > theseRegionsCollapsed.txt

  paste thisProm.txt theseRegionsCollapsed.txt >> RegionsCollapsed.txt

  #### COUNT PROCESSED GENES SO FAR
  counter=$((counter+1))
  if ! (($counter % 1000)); then
    echo "$counter promoters processed"
   # break
  fi
done < temp

awk -F\\t '{print $2}' RegionsCollapsed.txt | sed s/\,/\\n/g | sed '/^[[:space:]]*$/d'|sort | uniq | sed s/\_/\\t/g|sort -k1,1 -k2,2n > RegionsCollapsed.bed

##### CLEANUP

rm -rf thisProm.txt theseRegions.txt theseRegions.bed theseRegionsCollapsed.txt temp


#export iterator=1
#for arg in "${myArray[@]}"
#do
#  rm -rf temp$iterator
#  iterator=$((iterator+1))
#done
