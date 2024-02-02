dir=$1

module load bedtools/2.30.0

cd $dir

for i in pergene.pOtherEnds.{CTCF,enhancer,stringentEnhancer,H3K27me3,H3K4me3,promoter}.txt; do cut -f 2 $i | sed s/\,/\\n/g | sort | uniq |sed '/^[[:space:]]*$/d'|sed s/\_/\\t/g > ${i/.txt/.bed};done

CTCF=$(wc -l pergene.pOtherEnds.CTCF.bed|cut -d" " -f 1)

H3K27ac=$(wc -l pergene.pOtherEnds.enhancer.bed|cut -d" " -f 1)

Enhancer=$(wc -l pergene.pOtherEnds.stringentEnhancer.bed|cut -d" " -f 1)

promoter=$(wc -l pergene.pOtherEnds.promoter.bed|cut -d" " -f 1)

H3K27me3=$(wc -l pergene.pOtherEnds.H3K27me3.bed|cut -d" " -f 1)

H3K4me3=$(wc -l pergene.pOtherEnds.H3K4me3.bed|cut -d" " -f 1)

other=$(intersectBed  -v -a pergene.looped.bed -b pergene.pOtherEnds.enhancer.bed pergene.pOtherEnds.promoter.bed pergene.pOtherEnds.H3K4me3.bed pergene.pOtherEnds.CTCF.bed pergene.pOtherEnds.H3K27me3.bed|wc -l)

total=$(wc -l  pergene.looped.bed|cut -d" " -f 1)

CTCF_ratio=$(echo "scale=3;$CTCF/$total"|bc -l)
H3K27ac_ratio=$(echo "scale=3;$H3K27ac/$total"|bc -l)
Enhancer_ratio=$(echo "scale=3;$Enhancer/$total"|bc -l)
promoter_ratio=$(echo "scale=3;$promoter/$total"|bc -l)
H3K27me3_ratio=$(echo "scale=3;$H3K27me3/$total"|bc -l)
H3K4me3_ratio=$(echo "scale=3;$H3K4me3/$total"|bc -l)
other_ratio=$(echo "scale=3;$other/$total"|bc -l)

echo -e "total\t$total\nH3K27ac\t$H3K27ac\nH3K27ac_ratio\t$H3K27ac_ratio\nEnhancer\t$Enhancer\nEnhancer_ratio\t$Enhancer_ratio\nPromoter\t$promoter\nPromoter_ratio\t$promoter_ratio\nH3K4me3\t$H3K4me3\nH3K4me3_ratio\t$H3K4me3_ratio\nCTCF\t$CTCF\nCTCF_ratio\t$CTCF_ratio\nH3K27me3\t$H3K27me3\nH3K27me3_ratio\t$H3K27me3_ratio\nother\t$other\nother_ratio\t$other_ratio" >oEND.summary.txt






