dir=$1
TAD=$2

module load bedtools/2.30.0

cd $dir

inIN=$(grep -c 'inIN' summary.txt)
notIN=$(grep -c 'notIN' summary.txt)

hasINenh=$(grep -c 'hasINenh' summary.txt)
noINenh=$(grep -c 'noINenh' summary.txt)

hasLoopEnh=$(grep -c 'hasLoopEnh' summary.txt)
noLoopEnh=$(grep -c 'noLoopEnh' summary.txt)

let total=$inIN+$notIN

hasINenh_hasLoopEnh=$(grep -c 'hasINenh,hasLoopEnh' summary.txt)
hasINenh_noLoopEnh=$(grep -c 'hasINenh,noLoopEnh' summary.txt)

noINenh_hasLoopEnh=$(grep -c 'noINenh,hasLoopEnh' summary.txt)
noINenh_noLoopEnh=$(grep -c 'noINenh,noLoopEnh' summary.txt)


echo -e "IN+\t$inIN\nIN-\t$notIN\nINenh+\t$hasINenh\nINenh-\t$noINenh\nLoopEnh+\t$hasLoopEnh\nLoopEnh-\t$noLoopEnh\nLoopPromCoverage\t$LoopPromCoverage\ntotal\t$total" > gene_stat1.txt

echo -e "IN+Loop+\t$hasINenh_hasLoopEnh\nIN+Loop-\t$hasINenh_noLoopEnh\nIN-Loop+\t$noINenh_hasLoopEnh\nIN-Loop-\t$noINenh_noLoopEnh\ntotal\t$total" > gene_stat2.txt


enhInINCount=$(intersectBed -f 1 -u -wa -a theseEnhs.bed -b theseINs.bed | wc -l)
promInINCount=$(intersectBed -f 1 -u -wa -a theseProms.bed -b theseINs.bed | wc -l)
promUInINCount=$(intersectBed -f 1 -u -wa -a theseProms.uniqu.bed -b theseINs.bed | wc -l)
wholegeneInINCount=$(intersectBed -f 1 -u -wa -a theseGenes.bed -b theseINs.bed | wc -l)
genomeInIN=$(mergeBed -i theseINs.bed | awk -F\\t '{sum+=$3-$2} END {print sum}')


LoopEnhCoverage=$(intersectBed -u -wa -a theseEnhs.bed -b pairedbeds.bed | wc -l)
LoopPromCoverage=$(intersectBed -u -wa -a theseProms.bed -b pairedbeds.bed | wc -l)
LoopPromUCoverage=$(intersectBed -u -wa -a theseProms.uniqu.bed -b pairedbeds.bed | wc -l)
LoopWholegeneCoverage=$(intersectBed -u -wa -a theseGenes.bed -b pairedbeds.bed | wc -l)

enhInTAD=$(intersectBed -f 1 -u -a theseEnhs.bed -b $TAD|wc -l)
promInTAD=$(intersectBed -f 1 -u -a theseProms.bed -b $TAD|wc -l)
genomeInTAD=$(cut -f 1-3 $TAD |sort -k1,1 -k2,2n | mergeBed -i - | awk -F\\t '{sum+=$3-$2} END {print sum}')
genomeTotal=2725537669

enhByBoth=$(intersectBed -f 1 -wa -u -a theseEnhs.bed -b theseINs.bed|intersectBed -wa -u -a stdin -b pairedbeds.bed|wc -l)
promByBoth=$(intersectBed -f 1 -wa -u -a theseProms.bed -b theseINs.bed|intersectBed -wa -u -a stdin -b pairedbeds.bed|wc -l)

enhByLoopAndTAD=$(intersectBed -f 1 -wa -u -a theseEnhs.bed -b $TAD|intersectBed -wa -u -a stdin -b pairedbeds.bed|wc -l)
promByLoopAndTAD=$(intersectBed -f 1 -wa -u -a theseProms.bed -b $TAD|intersectBed -wa -u -a stdin -b pairedbeds.bed|wc -l)

enhByINAndTAD=$(intersectBed -f 1 -wa -u -a theseEnhs.bed -b $TAD|intersectBed -f 1 -wa -u -a stdin -b theseINs.bed|wc -l)
promByINAndTAD=$(intersectBed -f 1 -wa -u -a theseProms.bed -b $TAD|intersectBed -f 1 -wa -u -a stdin -b theseINs.bed|wc -l)

enhByThree=$(intersectBed -f 1 -wa -u -a theseEnhs.bed -b theseINs.bed|intersectBed -wa -u -a stdin -b pairedbeds.bed|intersectBed -f 1 -wa -u -a stdin -b $TAD|wc -l)
promByThree=$(intersectBed -f 1 -wa -u -a theseProms.bed -b theseINs.bed|intersectBed -wa -u -a stdin -b pairedbeds.bed|intersectBed -f 1 -wa -u -a stdin -b $TAD|wc -l)

let enhByEither=$enhInINCount+$LoopEnhCoverage-$enhByBoth
let promByEither=$promInINCount+$LoopPromCoverage-$promByBoth
let enhByEitherThree=$enhInINCount+$LoopEnhCoverage+$enhInTAD-$enhByBoth-$enhByLoopAndTAD-$enhByINAndTAD+$enhByThree
let promByEitherThree=$promInINCount+$LoopPromCoverage+$promInTAD-$promByBoth-$promByLoopAndTAD-$promByINAndTAD+$promByThree

totalEnh=$(wc -l theseEnhs.bed |cut -d" " -f 1)
totalProm=$(wc -l theseProms.bed |cut -d" " -f 1)
totalPromU=$(wc -l theseProms.uniqu.bed |cut -d" " -f 1)
totalIN=$(wc -l theseINs.bed |cut -d" " -f 1)
totalTAD=$(wc -l $TAD |cut -d" " -f 1)
totalLoop=$(wc -l pairedbeds.bed |cut -d" " -f 1)
totalLoop=$(echo $totalLoop/2|bc)
echo -e "enhInINCount\t$enhInINCount\npromInINCount\t$promInINCount\npromUInINCount\t$promUInINCount\nwholegeneInINCount\t$wholegeneInINCount\ngenomeInIN\t$genomeInIN\n\nLoopEnhCoverage\t$LoopEnhCoverage\nLoopPromCoverage\t$LoopPromCoverage\nLoopPromUCoverage\t$LoopPromUCoverage\nLoopWholegeneCoverage\t$LoopWholegeneCoverage\n\nenhInTAD\t$enhInTAD\npromInTAD\t$promInTAD\ngenomeInTAD\t$genomeInTAD\n\nenhByBoth\t$enhByBoth\npromByBoth\t$promByBoth\nenhByLoopAndTAD\t$enhByLoopAndTAD\npromByLoopAndTAD\t$promByLoopAndTAD\nenhByINAndTAD\t$enhByINAndTAD\npromByINAndTAD\t$promByINAndTAD\nenhByThree\t$enhByThree\npromByThree\t$promByThree\n\nenhByEither\t$enhByEither\npromByEither\t$promByEither\n\nenhByEitherThree\t$enhByEitherThree\npromByEitherThree\t$promByEitherThree\n\ntotalEnh\t$totalEnh\ntotalProm\t$totalProm\ntotalPromU\t$totalPromU\ntotalIN\t$totalIN\ntotalTAD\t$totalTAD\ntotalLoop\t$totalLoop\ntotalGenome\t$genomeTotal" > Loop_INs_TADs_stat.txt
#echo -e "enhInTAD\t$enhInTAD\npromInTAD\t$promInTAD\ngenomeInTAD\t$genomeInTAD\nenhTotal\t$enhTotal\npromTotal\t$promTotal\ngenomeTotal\t$genomeTotal" > TAD_stat.txt

cut -f 2 pergene.looped.txt | sed s/\,/\\n/g | sort | uniq |sed '/^[[:space:]]*$/d'|sed s/\_/\\t/g > pergene.looped.bed
cut -f 2 pergene.inIN.txt | sed s/\,/\\n/g | sort | uniq |sed '/^[[:space:]]*$/d'|sed s/\_/\\t/g > pergene.inIN.bed

#LoopEnhCoverage=$(intersectBed -u -wa -a theseEnhs.bed -b pairedbeds.bed | wc -l)
LoopEnh=$(intersectBed -u -wa -a theseEnhs.bed -b pergene.looped.bed | wc -l)
INenh=$(intersectBed -u -wa -a theseEnhs.bed -b pergene.inIN.bed | wc -l)

let LoopINenh=$totalEnh-$(wc -l EnhancersNotAssignableByLoops.bed|cut -d" " -f 1)
let LoopAndINenh=$LoopEnh+$INenh-$LoopINenh
let LoopUniq=$LoopEnh-$LoopAndINenh
let INUniq=$INenh-$LoopAndINenh
# total overlapped Enh (either overlap with Promoter or genebody)=Enh overlapped with genebody + Enh uniquely overlapped with promoters
let OverlapEnh=$( cut -f 2 NonLoopOverlappingGenes.txt |sed s/\,/\\n/g | sort | uniq |sed '/^[[:space:]]*$/d'|wc -l)+$(comm -13 <(cut -f 2 NonLoopOverlappingGenes.txt|sed s/\,/\\n/g | sort | uniq |sed '/^[[:space:]]*$/d') <(cut -f 2 NonLoopOverlappingProms.txt|sed s/\,/\\n/g|sort|uniq|sed '/^[[:space:]]*$/d') |wc -l)
Remainder=$(wc -l EnhancersNotAssignableByLoopsOverlap.bed |cut -d" " -f 1)

echo -e "LoopCoverage\t$LoopEnhCoverage\nLoopAssignable\t$LoopEnh\nIN\t$INenh\nLoop+IN+\t$LoopAndINenh\nLoopU\t$LoopUniq\nINU\t$INUniq\nLoop+IN\t$LoopINenh\nOverlap\t$OverlapEnh\nRemainder\t$Remainder\ntotal\t$totalEnh" >enh_stat.txt





