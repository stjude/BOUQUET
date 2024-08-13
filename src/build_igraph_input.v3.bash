
dir=$1
outdir=${2:-$dir}
ScriptDir=$3
CREmapping=$4
CREmapping_promSub=$5
EnSignal=${6:-'03172010_614ATAAXX_B5_mm9.sorted.bam'}
EnControl=${7:-'02202009_313D1AAXX_B3_mm9.sorted.bam'}
pergeneINFile=${8:-'pergene.inIN.txt'}
#/rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie/SRR1613251.sorted.bam
#/rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie/SRR1613249.sorted.bam

ID=$(basename $EnSignal)
module load bedtools/2.30.0
echo -e "Working Directory: $dir"
echo -e "Output Directory: $outdir"
echo -e "CREmapping: $CREmapping"
echo -e "CREmapping_promSub: $CREmapping_promSub"
echo -e "Signal Bam: $EnSignal"
echo -e "Control Bam: $EnControl"
echo -e "pregeneINFile: $pergeneINFile"

#exit
mkdir -p $outdir

cd $dir

#: <<'END'
# get nodes for igraph
#echo -e id"\t"Plabel"\t"Elabel"\t"name1"\t"name2 > $outdir/P.nodes.txt;sort -k1,1 -k2,2n -k4,4r theseProms.bed|awk '!seen[$1$2$3]++'|awk -F '[\t|]' -v OFS="\t" '{print $1"_"$2"_"$3,1,"P",$5,$4}' >> $outdir/P.nodes.txt
echo -e id"\t"Plabel"\t"Elabel"\t"name1"\t"name2 > $outdir/P.nodes.txt

join -j 1 -t $'\t' <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3"\t"$4}' theseProms.bed|sort -k1,1) <(awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3"\t"$5"\t"$6}' $CREmapping|sort -k1,1)|awk -F '[\t|]' -v OFS="\t"  '{print $1,$5,$4,$3,$2}' >> $outdir/P.nodes.txt

awk -F "\t" -v OFS="\t" '{print $1"_"$2"_"$3,0,1,"enhancer","enhancer"}' theseEnhs.bed > $outdir/E.nodes.txt

cat $outdir/P.nodes.txt $outdir/E.nodes.txt > $outdir/PE_nodes.txt

if [ ! -s $outdir/PE_nodes.Norm.$ID ]
    then
    echo -e "################# Quantify Enhancer signal for each individual region...###############\n"
    intersectBed -c -a <(tail -n +2 $outdir/PE_nodes.txt|awk -F ['\t'_] -v OFS='\t' '{print $1,$2,$3}') -b $EnSignal | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $outdir/PE_nodes.Signal.$ID

    echo -e "################# Quantify Input signal for each individual region...###############\n"
    intersectBed -c -a <(tail -n +2 $outdir/PE_nodes.txt|awk -F ['\t'_] -v OFS='\t' '{print $1,$2,$3}') -b $EnControl | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $outdir/PE_nodes.Input.$ID
    
    paste  $outdir/PE_nodes.Signal.$ID $outdir/PE_nodes.Input.$ID|awk -F '\t' -v OFS='\t' '{norm= ($2>=$4 ? $2-$4 : 0)} {print $1,norm}' > $outdir/PE_nodes.Norm.$ID
else
    echo -e "Detected per-CRE signal at $outdir/PE_nodes.Norm.$ID, skipping...\n"
fi
#END

#get (merged) promoter regions that need to be substracted 
cut -f4 $CREmapping_promSub|grep -v "._-1.-1"|tr "," "\n"|tr "_" "\t" > $outdir/promSub.bed

if [ ! -s $outdir/promSub.Norm.$ID ]
    then
    echo -e "################# Quantify Enhancer signal for each (merged)promoter region to be substracted...###############\n"
    intersectBed -c -a $outdir/promSub.bed  -b $EnSignal | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $outdir/promSub.Signal.$ID

    echo -e "################# Quantify Input signal for each (merged)promoter region to be substracted...###############\n"
    intersectBed -c -a $outdir/promSub.bed  -b $EnControl | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $outdir/promSub.Input.$ID
    
    paste  $outdir/promSub.Signal.$ID $outdir/promSub.Input.$ID|awk -F '\t' -v OFS='\t' '{norm= ($2>=$4 ? $2-$4 : 0)} {print $1,norm}' > $outdir/promSub.Norm.$ID
else
    echo -e "Detected per promSub signal at  $outdir/promSub.Norm.$ID , skipping...\n"
fi
python $ScriptDir/CRE_Signal.py -m $CREmapping_promSub -t $outdir/PE_nodes.Norm.$ID -p $outdir/promSub.Norm.$ID  > $outdir/PE_nodes.Both.$ID

head -1 $outdir/PE_nodes.txt |sed "s/$/\tSignal\tSignalSub/" > $outdir/PE_nodes.$ID.txt

join -t$'\t' -j 1 <(tail -n +2 $outdir/PE_nodes.txt |sort -k1,1) <(sort -k1,1 $outdir/PE_nodes.Both.$ID) >> $outdir/PE_nodes.$ID.txt

#END

# Get edges from loops for igraph
cat pairedbeds.bed |paste - -|awk -F "\t" -v OFS="\t" '{print $8,$4}' > $outdir/edges.xls
cat pairedbeds.bed |cut -f 1-3|sort -k1,1 -k2,2n|uniq > $outdir/loop_anchor.bed

intersectBed  -wa -wb -b theseEnhs.bed -a $outdir/loop_anchor.bed|awk -F "\t" -v OFS="\t" '{ print $1"_"$2"_"$3,$4"_"$5"_"$6}'|awk -F "\t" -v OFS="\t" '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x""a[x]}'|sort -t_ -k1,1 -k2,2n > $outdir/loop_anchor_enh.xls
intersectBed  -wa -wb -b theseProms.bed -a $outdir/loop_anchor.bed|awk -F "\t" -v OFS="\t" '{ print $1"_"$2"_"$3,$4"_"$5"_"$6";"$7}'|awk -F "\t" -v OFS="\t" '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x""a[x]}'|sort -t_ -k1,1 -k2,2n > $outdir/loop_anchor_pro.xls

cd $outdir

# correct a error of script directory

python $ScriptDir/IgraphParser.v1.0.py -l edges.xls -p loop_anchor_pro.xls -e loop_anchor_enh.xls

cd $dir

# Get edges from INs for igraph

# third field is type of edges: 0 for loop supported edges and 1 for INs supported edges

join -j1 -t $'\t' <(tail -n +2 $outdir/P.nodes.txt|awk -F "\t" -v OFS="\t" '{print $5"|"$4,$1}'|sort -k1,1) <(sort -k1,1 $pergeneINFile)|awk -F "\t" -v OFS="\t" '$3!="" {print $2,$3}'|tr "," "\t"|awk -F "\t" -v OFS="\t" '{for (i=2;i<=NF;i++) {print $1,$i,1,1} }' > $outdir/PE_edges_IN.txt

cat $outdir/PE_edges_loop.txt $outdir/PE_edges_IN.txt > $outdir/PE_edges.txt

#Clean Up

rm $outdir/PE_nodes.txt $outdir/PE_nodes.Signal.$ID $outdir/PE_nodes.Input.$ID $outdir/promSub.Signal.$ID $outdir/promSub.Input.$ID

#python ~/3DSE/script/IgraphParser.py -l edges.xls -p loop_anchor_pro.xls -e loop_anchor_enh.xls|grep  -v "#"|awk -F "\t" -v OFS="\t" '$1!=$2{print $0}'|awk -F '[\t_]' -v OFS="\t" '($3-$2)!=4000 && ($6-$5)==4000 {print $0}'|wc -l




