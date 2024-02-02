
module load bedtools/2.30.0
dir=$1
Enhancers=$2
Promoters=$3
WholeGenes=$4


mkdir -p $dir

cd $dir

promID=$(basename $Promoters)
enhID=$(basename $Enhancers)
wgID=$(basename $WholeGenes)
if [ ! -s $promID.mapped ]
then
#merge isoforms with same promoter, seprated by comma(,)

sort -k1,1 -k2,2n -k4,4 $Promoters|awk -F'\t' -v OFS="\t" '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$4 }END{ for(i in a) print i,a[i] }' > $promID.uniqueProm.txt


# merge overlaping enhancers/prmoters, seperated by semicolon(;). Also add two column indicating whether the resulted intervel is enhancer(col 5) and promoter(col 6)
cat <(cut -f 1-3 $Enhancers|awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$1"_"$2"_"$3}')  $promID.uniqueProm.txt|sort -k1,1 -k2,2n|mergeBed -i stdin -c 4 -o collapse -delim ";"|awk  -F "\t" -v OFS="\t"  -e '{ $5=($4 ~ /chr[[:alnum:]]+_[0-9]+/?1:0); $6=($4 ~ /(NM_)|(NR_)/?1:0);;print $0}' > $enhID.$promID.merged

# merge overlaping prmoters, seperated by semicolon(;). Also add two column indicating whether the resulted intervel is enhancer(col 5) and promoter(col 6)
cat $promID.uniqueProm.txt|sort -k1,1 -k2,2n|mergeBed -i stdin -c 4 -o collapse -delim ";" > $promID.merged

#Detect (mergered) promoter regions of each CRE ndoe, which will be removed when calculating exclusive enhancer signal loading per CRE
intersectBed -a $enhID.$promID.merged -b  $promID.merged  -wao|awk -F'\t' -v OFS="\t" '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$7"_"$8"_"$9 }END{ for(i in a) print i,a[i] }' > $enhID.$promID.merged.promSub

#generate new enhancer set
awk -F "\t" -v OFS="\t" '($5==1)&&($6==0){print $1,$2,$3,$4,$6}' $enhID.$promID.merged > $enhID.new
awk -F "\t" -v OFS="\t" '$5==1{print $1,$2,$3,$4,$6}' $enhID.$promID.merged > $enhID.new.inclusive
#generate new promoter set
#use sed to remove all enhancerID in col 4, in order to only keep the first NM/NR ID as representitive of each promoter
awk -F "\t" -v OFS="\t" '$6==1{print $1,$2,$3,$4}' $enhID.$promID.merged|sed -r 's/chr[[:alnum:]]*_[0-9]*_[0-9]*[,;]*//g'|sed -r 's/;$//'|awk -F '[\t,;]' -v OFS="\t" '{print $1,$2,$3,$4}' > $promID.new

#generate new whole gene set based on promoter set
join -1 4 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4 <(sort -k4,4 $WholeGenes) <(sort -k4,4 $promID.new)|sort -k1,1 -k2,2n > $wgID.new

# creat a mapping between merged gene and their reprensentitive gene
awk -F "\t" -v OFS="\t" '$6==1{print $1,$2,$3,$4}' $enhID.$promID.merged|sed -r 's/chr[[:alnum:]]*_[0-9]*_[0-9]*[,;]*//g'|sed -r 's/;$//'|cut -f 4|awk -F '[,;]' -v OFS="\t" '{ if (NF>1) {for(i=2; i<=NF; i++) a[$i]=$1 }} END{for (i in a) print i,a[i]}'|sort -k2,2 > $promID.mapped
fi
