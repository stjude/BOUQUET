#/bin/bash
#===============
# A pipeline which takes input of
# 1)promoter and gene annotation
# 2)Enahncer set
# 3)Loop set
# 4)Isoluted Neibourghoodloop set
# 5)Enhancer signal measurement(ChIP-Seq bam)

# to

# 1) assign (3D) enhancers to its target genes according to chromatin loops (Loops) and isulated neighborhoods(INs)
# 2) establish enhancer/promoter interaction networks
# 3) define gene regulatory landscape through community detection method, directly connected neighborhooh or daisy-chained neighborhood.
# 4) call 3DSE associated genes by ranking per-gene total enhancer signal

# Other utility includes:
# 1) assign linear enhancers to its target genes and further call Superenahcer associated genes by ranking the culmulative total enhancer signal per-gene (modified ROSE)
# 2) summarise the distribution of promoters and enhancers across different type of EP assignment evidence (eg. by loop, INs, overlap, etc...) *currently disabled*
# 3) compare SE genes defined by 3D and linear method *currently disabled*
# 4) annotate loop anchor with potential biological function (e.g., enhancer, promoter, isulator, etc...) *currently disabled*
#
# author: Jie Lu and Brian Abraham
# Abraham Lab
# Department of Computational Biology
# St. Jude Children’s Research Hospital
#===============

#===============
# sample execution command for this script:
# bash 3DSE.wraper.sh -o OutDir [-d ScriptDir] [ -p Promoters] [-I INs] [-l Loops] [-e Enhancers] [-s EnSignal] [-c EnControl]
#===============


#####Update for v5.0#############
#new switch -f
#new switch -F
# edge now have a new attribute: Loop supported only, IN supported only, Loop+IN supported.

#####Update for v7.0############
#fix the label of focal promoter and its direct neighbors, meaning forcelly including all direct neighours in the final community

#####Update for v8.0############
#v7 cause more problem than original method (v6). Add backfilling, which rescure any direct neighbours that are not included in the final community

#####Update for v10.0############
#added new statistics for efficiently pick positive genes by propotion of inter-genic enhancers and distal enhancers

#####Update for v11############
#added function of profiling promoter looped other ends in terms of enhancer,K4me3,CTCF,K27ac.

#####Update for v12############
###add functional characterization of promter interacting bins in terms of promoters in addition to peaks of H3K27ac, H3K4me3,H3K27me3,CTCF

#####Update for v16############
###pick INs based on promoter instead of full genebody

###Update for V17###########
#####set seed value a paramter, which used in R random process

###Update for V18#####
##1. In addtion to neighbors weighted LP with back filling (result under network/), now the pipeline also call community using Default LP with back filling (under netwrok.LP.default/)
##2. Added a new parameter (-O) to control whether to perform functional characterization of promoter-interacting regions (default: not perform this analysis).

###Update for V19####
#1. merge genes with unique promoter, seperated by comma(,)
#2. merge genes sharing the same merged promoter, seperated by semicolon. Only keep one gene representive.
#3. map each gene to its representative gene

###Update for V23###
#1.Performing E/P coverage analysis using strict Enhancer defination
#2.Merge genes into unique communities before drawing swoosh

###Update for V24###
#1.remove Enhancer singnal fully contained within promtoter regions
###Update for V25###
#add the ability to output a edgelist file with edge supporting types: 0=="Loop only"; 1=="IN only"; 2=="Both"
####Update for V26
# addtionally report community.E and community.order when forcing enhancer and promoter to be mutrually exclusive
##Update for V27
#perform TAD limited linear proximity based Enhancer Orphan rescure
#####Update for V28
#update CTCT bed file and Bam

#####Update for V29
#update 4 peak files: H3K27ac, H3K4me3, CTCF, H3K27me3
#use less stringent enhancer defination for characterizing promoter-interacted regions.
usage(){
    echo -e "usage : bash 3DSE.wraper.sh -o OutDir [-d ScriptDir] [-p Promoters] [-I INs] [-l Loops] [-e Enhancers] [-s EnSignal] [-c EnControl] [-g HighlightGenes] [-h]"
    echo -e "Use option -h for more information"

}


help(){
    usage;
    cat << EOF
    Options:
        -o  OutDir              user specified directory for all output of pipeline, recommond to reflet the unique combination of input files
        -d  ScriptDir           directory for all 3DSE scripts, default to be "/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script"
        -p  Promoters           list of gene promoters, in the format of Chr, Start, End, TranscriptID|GeneID, default to be 4Kb region around TSS of mm9
        -w  WholeGenes          list of whole gene body, in the format of Chr, Start, End, TranscriptID|GeneID, default to be whole gene body of mm9
        -I  INs                 Isolated Neibourghoods, in the format of Chr, Start, End, dault to be INs for mm9
        -l  Loops               Direct loops connecting two regulatory regions, in the format of bedpe. Default be to the ones defined from H3K27ac HiChIP from mESC.
        -e  Enhancers           Set of enhancers in the format of Chr, Start, End. Default be to the ones defined from H3K27ac ChIP-Seq from mESC.
        -s  EnSignal            ChiP-seq signal to quantify Enhancer signal instensity in number of reads. Default to be Bam file of MED1 ChIP sample from mESC.
        -c  EnControl           Control for ChiP-seq signal to quantify Enhancer signal instensity in number of reads. Default to be Bam file of MED1 Control(Input) sample from mESC.
        -g  HighlightGenes      a list of genes for which comprehenvise GRN figures and tracks are to be generated. Expecting format 'GeneSymbol1|RefSeqID1,GeneSymbol2|RefSeqID2...'. Please note the singal quote is needed because of the special charactar "|".
        -n  Network             a switch indicating whether to draw network figures (7) for highlighted genes. could be slow for genes with complex landscape.default to turn off.
        -O  OtherEnd            a switch indicating whether to perform functional characterization for promoter-interacting regions.
        -f  LO                  a switch indicating whether to use both loop and IN supported reguloary edges (default, 0) or only loop supported regulatory edges (1).
        -F  PromIN              a switch indicating whether to only assign E to P (default, 0) or assign both E and P to P (1) according to INs (default, 0).
        -v  H3K27acPeaks        H3K27ac peaks to define enahcners
        -x  H3K4me3Peaks        H3K4me3 peaks to define promoters
        -y  CTCFPeaks           CTCF peaks to define insulators
        -z  H3K27me3Peaks       H3K27me3 peaks to define PcG sites
        -E  ExpFile             geneSymbol-level expression file calculated based on RNA-Seq
        -S  Seed                Seed number of R random process. same seed value would  generate reproducible gene communities based on LP method.
        -t  TAD                 Topological Associated Domain used for rescue EO (Enhancer Orphan) after community based Enhancer assignment

EOF
exit
}

function fail {
    printf '%s\n' "$1" >&2  ## Send message to stderr. Exclude >&2 if you don't want it that way.
    exit "${2-1}"  ## Return a code specified by $2 or 1 by default.
}

if [ $# -lt 1 ]
then
    usage
    exit
fi

#==============
# default parameters
#=============

OutDir='OutDir'
ScriptDir='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script'
Promoters='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mm9.refseq.02_01_2017.4kbproms.bed.filtered'
WholeGenes='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mm9.refseq.02_01_2017.wholegenes.bed.filtered'
INs='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mES_INs_fromReview.nonProm.sorted.bed3'
Loops='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/HiChIP_mES_C3_UT_H3K27ac.Qvalue.0.01.H3K27ac_input_Peaks_1e-9.bed.mm9.refseq.02_01_2017.4kbproms.filtered.combined.5000.20000000.bedpe'
Enhancers='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K27ac_input_Peaks_1e-9.mm10.bed'
EnSignal='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/03172010_614ATAAXX_B5_mm9.sorted.bam'
EnControl='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/02202009_313D1AAXX_B3_mm9.sorted.bam'
HighlightGenes='Irf2bpl|NM_145836,Suz12|NM_199196,Hist1h4h|NM_153173'
Network=0
OtherEnd=0
LO=0
PromIN=0
H3K27acPeaks='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Res_H3K27ac-p9_kd-auto_peaks.bed'
H3K4me3Peaks='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Res_H3K4me3-p9_kd-auto_peaks.bed'
CTCFPeaks='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Res_CTCF_run-p9_kd-auto_peaks.bed'
H3K27me3Peaks='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Res_H3K27me3-W200-G400-FDR0.01-island.bed'
ExpFile='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mm10.expressed.gene.txt'
Seed=1
TAD='/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/MergedDixon2012_mm10.bed'
#==============

#==============
# pasre parameters from command-line input
#=============

while getopts ":o:d:p:w:I:l:e:s:c:g:nOfFv:x:y:z:E:S:t:h" opt;
    do
        case "$opt" in
            o) OutDir=$OPTARG;;
            d) ScriptDir=$OPTARG;;
            p) Promoters=$OPTARG;;
            w) WholeGenes=$OPTARG;;
            I) INs=$OPTARG;;
            l) Loops=$OPTARG;;
            e) Enhancers=$OPTARG;;
            s) EnSignal=$OPTARG;;
            c) EnControl=$OPTARG;;
            g) HighlightGenes=$OPTARG;;
            n) Network=1;;
            O) OtherEnd=1;;
            f) LO=1;;
            F) PromIN=1;;
            v) H3K27acPeaks=$OPTARG;;
            x) H3K4me3Peaks=$OPTARG;;
            y) CTCFPeaks=$OPTARG;;
            z) H3K27me3Peaks=$OPTARG;;
            E) ExpFile=$OPTARG;;
            S) Seed=$OPTARG;;
            t) TAD=$OPTARG;;
            h) help ;;
            \?) usage
                echo "error: unrecognized option -$OPTARG";
                exit 1
                ;;
            :)
                echo "Option -$OPTARG requires an argument." >&2
                usage
                exit 1
                ;;
        esac
    done


 if [[ -z $OutDir || -z $Promoters || -z $WholeGenes || -z $INs || -z $Loops || -z $Enhancers || -z $EnSignal || -z $EnControl || -z $ScriptDir || -z $HighlightGenes || -z $H3K27acPeaks || -z $H3K4me3Peaks || -z $CTCFPeaks || -z $H3K27me3Peaks || -z $ExpFile || -z $Seed || -z $TAD ]]; then
      usage
      exit
 fi

ID=$(basename $OutDir)

## Display input parameters

echo -e "################Input papameters#####################"
echo -e "OutDir: $OutDir"
echo -e "ID: $ID"
echo -e "ScriptDir: $ScriptDir"
echo -e "Promoters: $Promoters"
echo -e "WholeGenes: $WholeGenes"
echo -e "INs: $INs"
echo -e "Loops: $Loops"
echo -e "Enhancers: $Enhancers"
echo -e "Enhancer Signal: $EnSignal"
echo -e "Enhancer Signal Control: $EnControl"
echo -e "Highlight Genes: $HighlightGenes"
echo -e "flag for network figures: $Network"
echo -e "flag for functional analysis of promoter-interacting regions: $OtherEnd"
echo -e "flag for Loop Only edgess: $LO"
echo -e "flag for allowing P to P assignment: $PromIN"
echo -e "H3K27ac Peaks: $H3K27acPeaks"
echo -e "H3K4me3 Peaks: $H3K4me3Peaks"
echo -e "CTCF Peaks: $CTCFPeaks"
echo -e "H3K27me3 Peaks: $H3K27me3Peaks"
echo -e "expression file: $ExpFile"
echo -e "Seed value for R random process: $Seed"
echo -e "TADs: $TAD\n"

module load python/2.7.15-rhel7 bedtools/2.30.0 R/4.0.2
#merge enhancer and promoter, re-define enhancer, promoter and wholegene set
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Merge enhancer and promoter and re-define enhancer, promoter and wholegene set...###############\n"

promID=$(basename $Promoters)
enhID=$(basename $Enhancers)
wgID=$(basename $WholeGenes)


if [ ! -s $OutDir/$promID.mapped ]
then
    bash $ScriptDir/redefineEnhancerPromoter.v2.sh $OutDir $Enhancers $Promoters $WholeGenes
else
    echo -e "Detected merged promoters with mapping at $promID.mapped, skipping merging...\n"
fi

export Promoters_orig=$Promoters
export Enhancers_orig=$Enhancers
export WholeGenes_orig=$WholeGenes

export EnhancersInclusive="$OutDir/$enhID.new.inclusive"
export CREmapping="$OutDir/$enhID.$promID.merged"
export CREmapping_promSub="$OutDir/$enhID.$promID.merged.promSub"
export Enhancers="$OutDir/$enhID.new"
export Promoters="$OutDir/$promID.new"
export WholeGenes="$OutDir/$wgID.new"

echo -e "################Redefined CREmapping, Enhancers, Promoters and WholeGenes #####################"
echo -e "CREmapping: $CREmapping"
echo -e "CREmapping_promSub: $CREmapping_promSub"
echo -e "Enhancers: $Enhancers"
echo -e "EnhancersInclusive: $EnhancersInclusive"
echo -e "Promoters: $Promoters"
echo -e "WholeGenes: $WholeGenes\n"

#Asign Enhancers to target genes
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Asign Enhancers to target genes based on chromosomal 3D information (Loops and INs)...###############\n"

if [ ! -s $OutDir/RemainderClosest.txt ]
then
    bash $ScriptDir/assign_EP.v6.sh $Enhancers $Loops $INs $Promoters $WholeGenes $OutDir $ScriptDir $H3K27acPeaks $H3K4me3Peaks $CTCFPeaks $H3K27me3Peaks $Enhancers_orig $Promoters_orig
else
    echo -e "Detected 3D EP assignment results at $OutDir, skipping EP assignment...\n"
fi

cd $OutDir

############get summary statistics fro 3D SE results######
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Get summary statistics for Enhancer...###############\n"
if [ ! -s $OutDir/enh_stat.txt ]
then
    bash $ScriptDir/stats.v2.0.sh $PWD $TAD
else
    echo -e "Detected summary statistics at $OutDir, skipping summary statistics calculation...\n"
fi

############get summary statistics for promoter Other End distributions######
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Get summary statistics for promoter Other End...###############\n"
if [ $OtherEnd -eq 1 ]
then
#: << 'END'

if [ ! -s $OutDir/oEND.summary.txt ]
then
    bash $ScriptDir/oEND.v2.bash $PWD
else
    echo -e "Detected pOtherEnd distribution at $OutDir, skipping distribution calculation...\n"
fi


############Seperate other end to enhancer and non-enhancers######
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Get other end profile for enahcner and non-enhancer Other end seperately...###############\n"

##Process enhancer other end

if [ ! -s $OutDir/pergene.pOtherEnds.looped.enhancer/pergene.pOtherEnds.looped.enhancer.combined.xls ]
then

    mkdir -p pergene.pOtherEnds.looped.enhancer

    cp pergene.pOtherEnds.looped.enhancer.txt pergene.pOtherEnds.looped.enhancer

    cd pergene.pOtherEnds.looped.enhancer

    bash $ScriptDir/collapse_pergene_EP_lists_Jie.v2.0.sh pergene.pOtherEnds.looped.enhancer.txt

    for i in SRR1613251 SRR1613252 SRR299029 SRR1613255
    do
        intersectBed -c -a RegionsCollapsed.bed -b /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie_mm10/$i.bam | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $i.Signal
    done

    echo -e bin"\t"H3K27ac"\t"H3K4me3"\t"CTCF"\t"H3K27me3 > pergene.pOtherEnds.looped.enhancer.combined.xls;  paste SRR1613251.Signal SRR1613252.Signal SRR299029.Signal SRR1613255.Signal |cut -f 1,2,4,6,8 >> pergene.pOtherEnds.looped.enhancer.combined.xls
    #Clean Up
    rm SRR1613251.Signal SRR1613252.Signal SRR299029.Signal SRR1613255.Signal
else
    echo -e "Detected enhancer other end result at $OutDir, skipping...\n"
fi
## Process noEnhancer other end

cd $OutDir
if [ ! -s $OutDir/pergene.pOtherEnds.looped.noEnhancer/pergene.pOtherEnds.looped.noEnhancer.combined.xls ]
then

    mkdir -p pergene.pOtherEnds.looped.noEnhancer

    cp pergene.pOtherEnds.looped.noEnhancer.txt pergene.pOtherEnds.looped.noEnhancer

    cd pergene.pOtherEnds.looped.noEnhancer

    bash $ScriptDir/collapse_pergene_EP_lists_Jie.v2.0.sh pergene.pOtherEnds.looped.noEnhancer.txt

    for i in SRR1613251 SRR1613252 SRR299029 SRR1613255
    do
        intersectBed -c -a RegionsCollapsed.bed -b /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/archive/HiChIP/GSE62380/bowtie_mm10/$i.bam | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > $i.Signal
    done

    echo -e bin"\t"H3K27ac"\t"H3K4me3"\t"CTCF"\t"H3K27me3 > pergene.pOtherEnds.looped.noEnhancer.combined.xls;  paste SRR1613251.Signal SRR1613252.Signal SRR299029.Signal SRR1613255.Signal |cut -f 1,2,4,6,8 >> pergene.pOtherEnds.looped.noEnhancer.combined.xls

    #Clean Up
    rm SRR1613251.Signal SRR1613252.Signal SRR299029.Signal SRR1613255.Signal

else
    echo -e "Detected noEnhancer other end result at $OutDir, skipping...\n"
fi
#END
else
    echo -e "Functional analysis of promoter-interacting regions is turned off, skipping this process...\n"
fi
############get nodes and edges input for igraph######
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Get nodes and edges input for igraph ###############\n"
cd $OutDir
if [ ! -s $OutDir/igraph/PE_edges.txt ]
then
    if [ $PromIN -eq 1 ]
    then
        bash $ScriptDir/build_igraph_input.v3.bash $OutDir $OutDir/igraph $ScriptDir $CREmapping $CREmapping_promSub $EnSignal $EnControl pergene.inIN.PE.txt
    else
        bash $ScriptDir/build_igraph_input.v3.bash $OutDir $OutDir/igraph $ScriptDir $CREmapping $CREmapping_promSub $EnSignal $EnControl pergene.inIN.txt
    fi
else
    echo -e "Detected igraph input at $OutDir/igraph, skipping...\n"
fi

############ build E/P interaction network and infer Gene Regulatory Landscape######

echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Build E/P interaction network and infer Gene Regulatory Landscape  ###############\n"
BamID=$(basename $EnSignal)
if [ $LO -eq 1 ]
then
    EdgeFile=PE_edges_loop.txt
else
    EdgeFile=PE_edges.txt
fi
# Suppress "No such file or directory" error message for find

#if [[ !  -f "$(find $OutDir/igraph/network/network.plusIN.promoter.*.PE_nodes.$BamID.txt -type f -exec echo {} \; 2>/dev/null)" ]]
#then
if [ $Network -eq 1 ]
then
    Rscript $ScriptDir/GRN.preset.community.weightedLP.backfilling.v6.4.R $OutDir/igraph PE_nodes.$BamID.txt $EdgeFile $HighlightGenes $Seed
    Rscript $ScriptDir/GRN.preset.community.weightedLP.backfilling.v7.6.R $OutDir/igraph PE_nodes.$BamID.txt $EdgeFile $HighlightGenes $Seed
else
    Rscript $ScriptDir/GRN.preset.community.weightedLP.backfilling.v6.4.R $OutDir/igraph PE_nodes.$BamID.txt $EdgeFile
    Rscript $ScriptDir/GRN.preset.community.weightedLP.backfilling.v7.6.R $OutDir/igraph PE_nodes.$BamID.txt $EdgeFile
fi
#else
#    echo -e "Detected network result at $OutDir/igraph/network, skipping...\n"
#fi



########## Experimental module for generating ProteinPaint track, to be deployed ###############
: << 'END'
############ Generate ProteinPaint tracks for highlighted genes ######

echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Generate ProteinPaint tracks for highlighted genes  ###############\n"

if [[ ! "$(find $OutDir/igraph/network/*.community.enhancer.json -type f -exec echo {} \;)" ]]
then
    bash $ScriptDir/GenerateTrack.bash $OutDir/igraph/network $ScriptDir $Bedj $Enhancers $Promoters
else
    echo -e "Detected track for highlighted genes at $OutDir/igraph/network, skipping...\n"
fi
END

##hard-coded directory name, may change in the future##

network_file=$(find $OutDir/igraph/network/network.plusIN.promoter.*.PE_nodes.$BamID.txt -type f -exec echo {} \;)
network_file2=$(find $OutDir/igraph/network.LP.default/network.plusIN.promoter.*.PE_nodes.$BamID.txt -type f -exec echo {} \;)

outdir=$(dirname $network_file)
outdir2=$(dirname $network_file2)

#########Generate community statistics and linear assignment to rescue unassigned enhancers
#using stringent enhancer defination
bash $ScriptDir/community.stats.v5.sh $network_file $Enhancers $Promoters $OutDir/igraph/PE_nodes.$BamID.txt
bash $ScriptDir/community.stats.v7.sh $network_file2 $Enhancers $Promoters $OutDir/igraph/PE_nodes.$BamID.txt $TAD
#using inclusive enhancer defination
bash $ScriptDir/community.stats.v6.sh $network_file $EnhancersInclusive $Promoters $OutDir/igraph/PE_nodes.$BamID.txt
bash $ScriptDir/community.stats.v6.1.sh $network_file2 $EnhancersInclusive $Promoters $OutDir/igraph/PE_nodes.$BamID.txt

#######Update the final network table###########################
python $ScriptDir/RescueByLinearProximity.v4.py -n $network_file -a $outdir/linearAssigned.all.med1.txt -b  $outdir/linearAssigned.Genes.all.med1.txt -c  $outdir/linearAssigned.rest.all.med1.txt -e $outdir/linearAssigned.orig.all.med1.txt -E $ExpFile
python $ScriptDir/RescueByLinearProximity.v6.py -n $network_file2 -t $outdir2/RescuedByTAD.med1.txt -a $outdir2/linearAssigned.all.med1.txt -b  $outdir2/linearAssigned.Genes.all.med1.txt -c  $outdir2/linearAssigned.rest.all.med1.txt -e $outdir2/linearAssigned.orig.all.med1.txt -E $ExpFile

############ Call Super-enahcner associated genes based on Gene Regulatory landscape ######
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Call Super-enahcner associated genes based on E/P interaction network ###############\n"
if [ ! -s $network_file.update.uniqueCommunity.community.E.signal.AllGenes.table.txt ]
then
    bash $ScriptDir/GRN.callSuper.v3.bash $network_file.update $ScriptDir $HighlightGenes
    head -1 $network_file.update > $network_file.update.uniqueCommunitiy
    ##merge community with the same ID and same community.signal
    tail -n +2  $network_file.update|awk -F"\t" '!seen[$2, $32]++' >>  $network_file.update.uniqueCommunity
    tail -n +2  $network_file.update|awk -F"\t" -v OFS="\t" '{k=$2 FS $32; a[k] = a[k] == "" ? $1 : a[k] "," $1} END {for (k in a) print k,a[k]}' >  $network_file.update.uniqueCommunity.collapsed
    #mkdir -p $outdir/uniqueCommunity
    #cd $outdir/uniqueCommunity
    bash $ScriptDir/GRN.callSuper.v3.bash $network_file.update.uniqueCommunity $ScriptDir
else
    echo -e "Detected SuperGenes table $network_file.community.E.signal.SuperGenes.table.txt, skipping\n"
fi
##Call super-enhancer based on Default LP network with back-filling
if [ ! -s $network_file2.update.uniqueCommunity.community.E.signal.AllGenes.table.txt ]
then
    bash $ScriptDir/GRN.callSuper.v4.bash $network_file2.update $ScriptDir $HighlightGenes
    head -1 $network_file2.update > $network_file2.update.uniqueCommunity
    tail -n +2  $network_file2.update|awk -F"\t" '!seen[$2, $32]++' >>  $network_file2.update.uniqueCommunity
    tail -n +2  $network_file2.update|awk -F"\t" -v OFS="\t" '{k=$2 FS $32; a[k] = a[k] == "" ? $1 : a[k] "," $1} END {for (k in a) print k,a[k]}' >  $network_file2.update.uniqueCommunity.collapsed
    #mkdir -p $outdir2/uniqueCommunity
    #cd $outdir2/uniqueCommunity
    #need to set up HighlightGenes or a default mouse gene is used
    bash $ScriptDir/GRN.callSuper.v4.bash $network_file2.update.uniqueCommunity $ScriptDir
else
    echo -e "Detected SuperGenes table $network_file2.community.E.signal.SuperGenes.table.txt, skipping\n"
fi

########the old 3DSE method without establish E/P interaction network##############
######### Collapse Enhancers assign by different ways to single gene
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Collapse Enhancers assigned by different ways to single gene...###############\n"
if [ ! -s RegionsCollapsed.txt ]
then
    bash $ScriptDir/collapse_pergene_EP_lists_Jie.v2.0.sh pergene.looped.txt pergene.inIN.txt NonLoopOverlappingProms.txt NonLoopOverlappingGenes.txt RemainderClosest.txt
else
    echo -e "Detected pergnee collapsed EP assignment result at RegionsCollapsed.txt, skipping EP assignment collapse...\n"
fi
# Remove interchromosomal Enhancer-Target assignment (currently not functional, since interchromosomal loops are not allowed at the first place)
echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Remove interchromosomal Enhancer-Target assignment...###############\n"
if [ ! -s RegionsCollapsed.cis.txt ]
then
    bash $ScriptDir/cis_filter.sh RegionsCollapsed.txt $Promoters RegionsCollapsed.cis.txt $ScriptDir
else
    echo -e "Detected intra-chromosomal EP assignment only  result at RegionsCollapsed.cis.txt, skipping interchromosomal EP  assignment removal...\n"
fi
# Quantify signal for each individual enhancer region. This process resquests a fairly large memory when intersecting with bam file, 30G is working for ~900MB bam file, be careful!!, added by Jie, 11/03/2020

echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Quantify ChIP signal for each individual enhancer region...###############\n"
if [ ! -s RegionsCollapsed.bed.Input ]
then
    intersectBed -c -a RegionsCollapsed.bed -b $EnSignal | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > RegionsCollapsed.bed.Signal

    echo -e  "--------------------------------------------"
    echo $(date)
    echo -e "################# Quantify Input signal for each individual enhancer region...###############\n"
    intersectBed -c -a RegionsCollapsed.bed -b $EnControl | awk -F\\t '{print $1 "_" $2 "_" $3 "\t" $4}' > RegionsCollapsed.bed.Input
else
    echo -e "Detected per-enhancer signal at RegionsCollapsed.bed.Signal and RegionsCollapsed.bed.Input, skipping...\n"
fi
# Quantify the total Enhancer signal for single transcript/gene.

echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Quantify the total Enhancer signal for single transcript/gene...###############\n"

if [ ! -s pergene.signal.txt ]
then
    perl $ScriptDir/synthesize_signal.pl pergene.signal.txt RegionsCollapsed.txt RegionsCollapsed.bed.Signal RegionsCollapsed.bed.Input
else
    echo -e "Detected pergene enhancer signal at pergene.signal.txt, skipping...\n"
fi

echo -e  "--------------------------------------------"
echo $(date)
echo -e "################# Call 3D SE...###############\n"


if [ ! -s $ID.AllGenes.table.txt ]
then
    R --save $ID pergene.signal.txt $HighlightGenes < $ScriptDir/PriMROSE_callSuper.flexibleID.R
else
    echo -e "Detected 3DSE results at $ID.AllGenes.table.txt, skipping...\n"
fi
