network_file=$1
enhSet=$2
promSet=$3
EnhSignlLoading=$4

module load bedtools/2.30.0

dir=$(dirname $network_file)/EnhancersInclusive
mkdir -p $dir
cd $dir
enhAssigned=$(comm -12 <(tail -n +2 $network_file|cut -f 1,2,3,4,7,25,26,27,33|awk -F "\t" -v OFS="\t" '$6>0 {print $9}'|tr ',' "\n"|sort|uniq) <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3}' $enhSet|sort -k1,1)|wc -l)
PromAssigned=$(tail -n +2 $network_file|cut -f 1,2,3,4,7,25,26,27,33|awk -F "\t" -v OFS="\t" '$6>0 {print $9}'|wc -l)
PromTotal=$(tail -n +2 $network_file|wc -l)
enhTotal=$(cat $enhSet|wc -l)

awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseProms.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a $enhSet -b - -d|cut -f 1,2,3,4,9,10|sort -k6,6n > linearAssigned.orig.all.txt
join -j1 -t $'\t' <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}' linearAssigned.orig.all.txt|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1,$6}' $EnhSignlLoading|sort -k1,1) > linearAssigned.orig.all.med1.txt

comm -12 <(tail -n +2 $network_file|cut -f 1,2,3,4,7,25,26,27,33|awk -F "\t" -v OFS="\t" '$6>0 {print $9}'|tr ',' "\n"|sort|uniq) <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3}' $enhSet|sort -k1,1)|sed '/^[[:space:]]*$/d'|sed s/\_/\\t/g > enhAssigned.bed
#join -1 4 -2 1 -t $'\t' <(sort -t $'\t' -k4,4 $promSet) <(tail -n +2 $network_file|cut -f 1,2,3,4,5,23,24,25,31|awk -F "\t" -v OFS="\t" '$6>0 {print $1}'|sort) > promAssigned.bed
join -j 1  -t $'\t' <(awk -F '[\t|]' -v OFS="\t" '{print $5"|"$4,$1,$2,$3}' $promSet|sort  -t $'\t' -k1,1) <(tail -n +2 $network_file|cut -f 1,2,3,4,7,25,26,27,33|awk -F "\t" -v OFS="\t" '$6>0 {print $1}'|sort)|awk -F "\t" -v OFS="\t" '{print $2,$3,$4,$1}' > promAssigned.bed

intersectBed -v -a $enhSet -b enhAssigned.bed > enhNotAssigned.bed

###rescue by promoter
awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseProms.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a enhNotAssigned.bed -b - -d|cut -f 1,2,3,4,9,10|sort -k6,6n > linearAssigned.all.txt
awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseProms.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a enhNotAssigned.bed -b - -d|cut -f 1,2,3,4,9,10|sort -k4,4 -u|sort -k6,6n > linearAssigned.txt
join -j1 -t $'\t' <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}' linearAssigned.all.txt|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1,$6}' $EnhSignlLoading|sort -k1,1) > linearAssigned.all.med1.txt
join -j1 -t $'\t' <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}' linearAssigned.all.txt|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1,$6}' $EnhSignlLoading|sort -k1,1)|sort -k1,1 -u > linearAssigned.med1.txt
awk -F "\t" -v OFS="\t" '$6 <= 5000 {print $0}' linearAssigned.txt > linearAssigned.0-5k.txt
cut -f 1-3 linearAssigned.0-5k.txt|sort -k1,1 -k2,2n > enhAssigned.promoter5k.bed
awk -F "\t" -v OFS="\t" '$6 <= 5000 {print $0}' linearAssigned.all.txt |cut -f 5|sort|uniq > linearAssigned.0-5k.promoters.txt
awk -F "\t" -v OFS="\t" '$6 > 5000 {print $0}' linearAssigned.txt > linearAssigned.5k-inf.txt
cut -f 1-4 linearAssigned.5k-inf.txt|sort -k1,1 -k2,2n > enhNotAssigned.distal.promoter5k.bed
cut -f 1-3 linearAssigned.5k-inf.txt|sort -k1,1 -k2,2n > enhNotAssigned.afterPromoter5kRescue.bed
cat <(cut -f 4 promAssigned.bed)  linearAssigned.0-5k.promoters.txt |sort|uniq > promAssignedBycommPlusLinerProm.txt
###rescue by genebody
awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseGenes.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a enhNotAssigned.distal.promoter5k.bed -b - -d|cut -f 1,2,3,4,8,9|sort -k4,4 -u|sort -k6,6n > linearAssigned.Genes.txt
awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseGenes.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a enhNotAssigned.distal.promoter5k.bed -b - -d|cut -f 1,2,3,4,8,9|sort -k6,6n > linearAssigned.Genes.all.txt
join -j1 -t $'\t' <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}' linearAssigned.Genes.all.txt|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1,$6}' $EnhSignlLoading|sort -k1,1) > linearAssigned.Genes.all.med1.txt
awk -F "\t" -v OFS="\t" '$6 <= 5000 {print $0}' linearAssigned.Genes.txt > linearAssigned.Genes.0-5k.txt
awk -F "\t" -v OFS="\t" '$6 <= 5000 {print $0}' linearAssigned.Genes.all.txt |cut -f 5|sort|uniq > linearAssigned.Genes.0-5k.promoters.txt
awk -F "\t" -v OFS="\t" '$6 > 5000 {print $0}' linearAssigned.Genes.txt > linearAssigned.Genes.5k-inf.txt
cut -f 1-4 linearAssigned.Genes.5k-inf.txt|sort -k1,1 -k2,2n > enhNotAssigned.distal.Genes5k.bed

join -j1 -t $'\t' <(awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseProms.uniqu.bed|sort  -t $'\t' -k1,1 -k2,2n|closestBed -t all -a enhNotAssigned.distal.Genes5k.bed -b - -d|cut -f 1,2,3,4,8,9|sort -k6,6n|awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}'|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1,$6}' $EnhSignlLoading|sort -k1,1) > linearAssigned.rest.all.med1.txt

#join -j 1  -t $'\t' <(awk -F '[\t|]' -v OFS="\t" '{print $5"|"$4,$1,$2,$3}' $promSet|sort  -t $'\t' -k1,1) <(tail -n +2 $network_file|cut -f 1,2,3,4,5,23,24,25,31|awk -F "\t" -v OFS="\t" '$6>0 {print $1}'|sort)|awk -F "\t" -v OFS="\t" '{print $2,$3,$4,$1}' > promAssigned.bed



PromRescued_5k=$(comm -13 <(cut -f 4 promAssigned.bed|sort) linearAssigned.0-5k.promoters.txt |wc -l)
EnhRescued_5K=$(cat linearAssigned.0-5k.txt |wc -l)

PromRescued_Genebody=$(comm -12 <(awk -F '[\t|]' -v OFS="\t" '{print $1,$2,$3,$5"|"$4}' ../../../theseProms.uniqu.bed|cut -f 4|sort -t $'\t' -k1,1)  <(comm -13 promAssignedBycommPlusLinerProm.txt linearAssigned.Genes.0-5k.promoters.txt) |wc -l)
EnhRescued_Genebody=$(cat linearAssigned.Genes.0-5k.txt |wc -l)

echo -e "enhAssigned\t$enhAssigned\nPromAssigned\t$PromAssigned\nPromTotal\t$PromTotal\nenhTotal\t$enhTotal\nEnhRescued_5K\t$EnhRescued_5K\nPromRescued_5k\t$PromRescued_5k\nEnhRescued_Genebody\t$EnhRescued_Genebody\nPromRescued_Genebody\t$PromRescued_Genebody" > community_stats.txt
