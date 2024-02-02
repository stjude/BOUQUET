#####compared to v2.0, instead of all looped bins, only consider pre-defined enhancers that overlap with a looped pOotherEnds bin. this greatly reduced the # of candidate enhancers considered downstream,  added by Jie 02/09/2021
###compared to v3, add functional characterization based on promoters instead of peaks of H3K27ac, H3K4me3,H3K27me3,CTCF
###copared to v4, use promoter regions rather than whole gene body when to pick the smallest INs for a given gene.
export enhancerSet=$1
export loopConnectivitySet=$2    # EXPECT CHR1\tSTART1\tEND1\tCHR2\tSTART2\tEND2
export insulatedNeighborhoods=$3 # EXPECT CHR\tSTART\tEND
export promoters=$4              # EXPECT CHR\tSTART-2kb\tSTART+2kb\tGENEID|GENENAME
export wholegenes=$5             # EXPECT CHR\tSTART\tEND\tGENEID|GENENAME
export output=$6
export ScriptDir=$7
export H3K27acPeaks=$8
export H3K4me3Peaks=$9
export CTCFPeaks=${10}
export H3K27me3Peaks=${11}
export Enhancers_orig=${12}
export Promoters_orig=${13}


echo -e "\n\n"
echo -e "script:\t$0\nenahncerSet:\t$1\nloopConnectivitySet:\t$2\ninsulatedNeighborhoods:\t$3\npromoters:\t$4\nwholegenes:\t$5\noutput:\t$6\nScriptDir:\t$7\nH3K27acPeaks:\t$8\nH3K4me3Peaks:$9\nCTCFPeaks:\t${10}\nH3K27me3Peaks:\t${11}"
echo -e "\n\n"
echo -e "=========ENHANCER SET========="
ls $enhancerSet
echo -e "=====LOOPS FOR CONNECTING====="
ls $loopConnectivitySet
echo -e "====INSULATED NEIGHBORHOODS==="
ls $insulatedNeighborhoods
echo -e "===========PROMOTERS=========="
ls $promoters
echo -e "\n"

module load bedtools/2.30.0

mkdir -p $output
cd $output
#### PREPROCESS NEIGHBORHOODS ####
sort -k1,1 -k2,2n $insulatedNeighborhoods > $output/theseINs.bed

awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$3-$2}' $insulatedNeighborhoods|sort -k4,4n|cut -f 1-3 > $output/theseINs.sorted.bed

#keep a version of INs sorted by increasing size, which is used to choose the smallest IN for promoter

export insulatedNeighborhoods="$output/theseINs.sorted.bed"


#### PREPROCESS PROMOTERS ####
awk -F\\t '{print $4}' $promoters > $output/prom.txt
sort -k1,1 -k2,2n $promoters > $output/theseProms.bed
##get unique promoters
sort -k1,1 -k2,2n -k4,4r $output/theseProms.bed|awk '!seen[$1$2$3]++'|awk -F '\t' -v OFS="\t" '{print $0}' > theseProms.uniqu.bed
export promotersU="$output/theseProms.uniqu.bed"
export promoters="$output/theseProms.bed"

#### PREPROCESS WHOLEGENE ####
sort -k1,1 -k2,2n $wholegenes > $output/theseGenes.bed
sort -k1,1 -k2,2n -k4,4r $output/theseGenes.bed|awk '!seen[$1$2$3]++'|awk -F '\t' -v OFS="\t" '{print $0}' > theseGenes.uniqu.bed
export genebodys="$output/theseGenes.bed"

#### PREPROCESS ENHANCERS ####
## make sure only the first 3 columns (chr, start, end) are kept, modified by Jie 10/26/2020
## only keep enahcners that are not overlap with promoters, modified by Jie 07/27/2021
#intersectBed -a $enhancerSet -b $promoters -v|sort -k1,1 -k2,2n|cut -f 1-3  > $output/theseEnhs.bed
##keep all full enhancers,noting excluding the ones overlaping with promoters,  modified by Jie 05/31/2022
sort -k1,1 -k2,2n $enhancerSet |cut -f 1-3  > $output/theseEnhs.bed
export enhancerSet="$output/theseEnhs.bed"

#### MAKE USEFUL VISUALIZATIONS
awk -F\\t '{print $1 ":" $2 "-" $3 "\t" $4 ":" $5 "-" $6 "\t1"}' $loopConnectivitySet > $output/connLoops.washu

#### GET USEFUL NUMBERS
export enhInINCount=$(intersectBed -u -wa -a $enhancerSet -b $insulatedNeighborhoods | wc -l)
export promInINCount=$(intersectBed -u -wa -a $promoters -b $insulatedNeighborhoods | wc -l)
export promUInINCount=$(intersectBed -u -wa -a $promotersU -b $insulatedNeighborhoods | wc -l)
export wholegeneInINCount=$(intersectBed -f 1 -u -wa -a $wholegenes -b $insulatedNeighborhoods | wc -l)
export enhCount=$(wc -l $enhancerSet | awk -F' ' '{print $1}')
export promCount=$(wc -l ${output}/theseProms.bed | awk -F' ' '{print $1}')
export promCountUniq=$(cut -f 1-3 $output/theseProms.bed |uniq|wc -l| awk -F' ' '{print $1}')
export INcount=$(wc -l $insulatedNeighborhoods | awk -F' ' '{print $1}')
export genomeInIN=$(mergeBed -i $output/theseINs.bed | awk -F\\t '{sum+=$3-$2} END {print sum}')

echo -e "PROM COUNT:\t$promCount"
echo -e "UNIQUE PROM COUNT:\t$promCountUniq"
echo -e "ENHANCER COUNT:\t$enhCount"
echo -e "INSULATED NEIGHBORHOOD COUNT:\t$INcount"
echo -e "GENOME IN NEIGHBORHOODS:\t$genomeInIN"
echo -e "ENHANCER IN NEIGHBORHOODS\t$enhInINCount"
echo -e "PROM IN NEIGHBORHOODS:\t$promInINCount"
echo -e "Unique PROM IN NEIGHBORHOODS:\t$promUInINCount"
echo -e "WHOLEGENE IN NEIGHBORHOODS:\t$wholegeneInINCount"
echo -e "\n"

#### FIND PROMOTERS THAT HAVE CONNECTING LOOPS
awk -F\\t '{print $1 "_" $2 "_" $3 "\n" $4 "_" $5 "_" $6}' $loopConnectivitySet | sort | uniq | sed s/\_/\\t/g > $output/connLoopAnchors.bed3
intersectBed -loj -a $output/connLoopAnchors.bed3 -b $promoters > $output/connLoopPromAnchors.bed3

#########
##for loop anchor without overlapping promoters, there are 8 fileds instead of 7, because there is a empty string between -1($6) and .($8).
## corrected by Jie 10/20/2020
##########

awk -F\\t '{if($7 != "") print $7  "\t" $1 "_" $2 "_" $3}' $output/connLoopPromAnchors.bed3 > $output/connLoopPromAnchors.txt
#awk -F\\t '{if($7 != ".") print $7  "\t" $1 "_" $2 "_" $3}' $output/connLoopPromAnchors.bed3 > $output/connLoopPromAnchors.txt
awk -F\\t '{print $1}' $output/connLoopPromAnchors.txt | sort | uniq > $output/promsWithConnLoops


#### BUILD ANCHOR-ANCHOR FILE
awk -F\\t '{print $1 "\t" $2 "\t" $3 "\t" $4 "_" $5 "_" $6 "\n" $4 "\t" $5 "\t" $6 "\t" $1 "_" $2 "_" $3}' $loopConnectivitySet > $output/pairedbeds.bed


#### SET UP OUTPUT FILES
for file in  $output/pergene.looped.txt $output/pergene.inIN.txt $output/pergene.INs.txt $output/pergene.inIN.PE.txt $output/log $output/summary.txt $output/pergene.IN.txt $output/pergene.pOtherEnds.looped.enhancer.txt $output/pergene.pOtherEnds.looped.noEnhancer.txt $output/pergene.pOtherEnds.enhancer.txt $output/pergene.pOtherEnds.H3K4me3.txt $output/pergene.pOtherEnds.CTCF.txt $output/pergene.pOtherEnds.H3K27me3.txt $output/pergene.pOtherEnds.promoter.txt
#for file in  $output/pergene.looped.txt $output/pergene.inIN.txt $output/log $output/summary.txt $output/pergene.IN.txt $output/pergene.pOtherEnds.enhancer.txt $output/pergene.pOtherEnds.H3K4me3.txt $output/pergene.pOtherEnds.CTCF.txt $output/pergene.pOtherEnds.H3K27me3.txt $output/pergene.pOtherEnds.looped.enhancer.txt $output/pergene.pOtherEnds.looped.noEnhancer.txt

do
    [ ! -e $file ] || rm $file
    touch $file
done


export counter=1

#### FOR EACH PROMOTER, FIND ENHANCERS ASSIGNABLE MULTIPLE WAYS
while read thisProm
do
  echo -e "===$thisProm===" >> $output/log

### ASSIGN REGIONS DIRECTLY LOOPED
  grep -w $thisProm $promoters | awk -F\\t '{print $1 "\t" $2 "\t" $3}' > $output/thisProm.bed
  grep -w $thisProm $genebodys | awk -F\\t '{print $1 "\t" $2 "\t" $3}' > $output/thisGene.bed
  # get uniq loop anchors, added by jie 01/19/2021

  intersectBed -b $output/thisProm.bed -a $output/pairedbeds.bed | awk -F\\t '{print $4}' | sed s/\_/\\t/gi|sort -k1,1 -k2,2n|uniq > $output/pOtherEnds.looped.all.bed
  # a alternatave version, instead of all looped bins, only consider pre-defined enhancers that overlap with a looped pOotherEnds bin, added by Jie 02/09/2021
  intersectBed -wa -u -b $output/pOtherEnds.looped.all.bed -a $enhancerSet|sort -k1,1 -k2,2n  > $output/pOtherEnds.looped.bed

  #divide the looped pOtherEnds bins into the ones overlapping with enhancer and the ones not.
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $enhancerSet  > $output/pOtherEnds.looped.enhancer.bed

  intersectBed -v -a $output/pOtherEnds.looped.all.bed -b $enhancerSet  > $output/pOtherEnds.looped.noEnhancer.bed

  #intersectBed -v -a $output/pOtherEnds.looped.bed -b $output/thisProm.bed | sort -k1,1 -k2,2n > $output/pOtherEnds.nonProm.looped.bed
  ##set lable a(:a), accumulate all lines into pattern space(N),branching to a (ba) until the end of file ($!); replace newline to comma globally annoatated by Jie 10/26/2020  #####
  sed ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.looped.all.bed > $output/pOtherEnds.looped.all.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.looped.all.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.looped.all.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.looped.txt

  sed   ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.looped.enhancer.bed > $output/pOtherEnds.looped.enhancer.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.looped.enhancer.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.looped.enhancer.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.looped.enhancer.txt

  sed   ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.looped.noEnhancer.bed > $output/pOtherEnds.looped.noEnhancer.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.looped.noEnhancer.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.looped.noEnhancer.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.looped.noEnhancer.txt

  ### ASSIGN REGIONS IN SAME IN
# INs are sorted by size, therefore "head -n 1" is to make sure to choose the smallest INs among multiple ones
  #intersectBed -F 1 -wa -u -a $insulatedNeighborhoods -b $output/thisGene.bed  > $output/theseINs.bed3
# instead of full genebody, use promoter regions to select INs for a given gene
  intersectBed -F 1 -wa -u -a $insulatedNeighborhoods -b $output/thisProm.bed  > $output/theseINs.bed3
  head -n 1 $output/theseINs.bed3 > $output/thisIN.bed3
# mearge all nested INs for a promoter,calculating the largest possible encompassing IN region
  sort -k1,1 -k2,2n $output/theseINs.bed3 > $output/theseINs.sorted.bed3
  mergeBed -i $output/theseINs.sorted.bed3 >$output/theseINs.merged.bed3

  intersectBed -u -a $enhancerSet -b $output/thisIN.bed3 > $output/pOtherEnds.inIN.bed
  intersectBed -v -a $output/pOtherEnds.inIN.bed -b $output/thisProm.bed | sort -k1,1 -k2,2n > $output/pOtherEnds.nonProm.inIN.bed

  sed  ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.inIN.bed >$output/pOtherEnds.inIN.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.inIN.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.inIN.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.inIN.txt
  #record IN for each promoter
  paste $output/line $output/theseINs.bed3 >> $output/pergene.INs.txt

  # also assign promoters within the INs to the target promoters
  intersectBed -u -a $promotersU -b $output/thisIN.bed3|cut -f 1-3 > $output/pOtherEnds.prom.inIN.bed
  cat  $output/pOtherEnds.inIN.bed $output/pOtherEnds.prom.inIN.bed >  $output/pOtherEnds.inIN.both.bed

  sed  ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.inIN.both.bed >$output/pOtherEnds.inIN.both.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.inIN.both.bed2
  paste $output/line $output/pOtherEnds.inIN.both.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.inIN.PE.txt

#: <<'END'
  #echo $thisProm > $output/line
  #paste $output/line $output/thisIN.bed3 >> $output/pergene.IN.txt
  ## COUNT of distribution of pOtherEnds in enhancers(K3k27 peaks), promoters(k4me3 peaks), insulators (CTCF peaks) and PcG sites (K27me3)
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $H3K27acPeaks > $output/pOtherEnds.enhancer.bed
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $H3K4me3Peaks > $output/pOtherEnds.H3K4me3.bed
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $CTCFPeaks > $output/pOtherEnds.CTCF.bed
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $H3K27me3Peaks > $output/pOtherEnds.H3K27me3.bed
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.bed -b $promoters > $output/pOtherEnds.promoter.bed
  # get promoter interacted bins that not overlapping with other promoters
  intersectBed -a $output/pOtherEnds.looped.all.bed -b $Promoters_orig -v > $output/pOtherEnds.looped.all.nonProm.bed
  # get the non-promoters bins that overlap with H3K27ac peaks (enhancers).
  intersectBed -wa -u -a $output/pOtherEnds.looped.all.nonProm.bed -b $Enhancers_orig > $output/pOtherEnds.stringentEnhancer.bed

  sed   ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.enhancer.bed > $output/pOtherEnds.enhancer.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.enhancer.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.enhancer.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.enhancer.txt

  sed   ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.stringentEnhancer.bed > $output/pOtherEnds.stringentEnhancer.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.stringentEnhancer.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.stringentEnhancer.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.stringentEnhancer.txt

  sed  ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.H3K4me3.bed > $output/pOtherEnds.H3K4me3.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.H3K4me3.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.H3K4me3.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.H3K4me3.txt

  sed   ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.CTCF.bed > $output/pOtherEnds.CTCF.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.CTCF.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.CTCF.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.CTCF.txt

  sed  ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.H3K27me3.bed > $output/pOtherEnds.H3K27me3.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.H3K27me3.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.H3K27me3.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.H3K27me3.txt

  sed  ':a;N;$!ba;s/\n/,/g' $output/pOtherEnds.promoter.bed > $output/pOtherEnds.promoter.bed2
  sed -i s/\\t/\_/g $output/pOtherEnds.promoter.bed2
  echo $thisProm > $output/line
  paste $output/line $output/pOtherEnds.promoter.bed2 > $output/line2
  cat $output/line2 >> $output/pergene.pOtherEnds.promoter.txt

  if [ -s $output/pOtherEnds.enhancer.bed ];then
    pOtherEnds_enhancer=$(wc -l $output/pOtherEnds.enhancer.bed|cut -d " " -f 1)
  else
    pOtherEnds_enhancer=0
  fi

  if [ -s $output/pOtherEnds.stringentEnhancer.bed ];then
    pOtherEnds_StringentEnhancer=$(wc -l $output/pOtherEnds.stringentEnhancer.bed|cut -d " " -f 1)
  else
    pOtherEnds_StringentEnhancer=0
  fi

  if [ -s $output/pOtherEnds.H3K4me3.bed ];then
    pOtherEnds_H3K4me3=$(wc -l $output/pOtherEnds.H3K4me3.bed|cut -d " " -f 1)
  else
    pOtherEnds_H3K4me3=0
  fi

  if [ -s $output/pOtherEnds.CTCF.bed ];then
    pOtherEnds_CTCF=$(wc -l $output/pOtherEnds.CTCF.bed|cut -d " " -f 1)
  else
    pOtherEnds_CTCF=0
  fi

  if [ -s $output/pOtherEnds.H3K27me3.bed ];then
    pOtherEnds_H3K27me3=$(wc -l $output/pOtherEnds.H3K27me3.bed|cut -d " " -f 1)
  else
    pOtherEnds_H3K27me3=0
  fi

  if [ -s $output/pOtherEnds.promoter.bed ];then
    pOtherEnds_promoter=$(wc -l $output/pOtherEnds.promoter.bed|cut -d " " -f 1)
  else
    pOtherEnds_promoter=0
  fi
#END
  ## COUNT LOOPED REGULATARY ELEMENTS (BINS) THAT ESCAPE INs

  # count the number of looped enhancers that in the same IN
  intersectBed -c  -F 1 -b $output/pOtherEnds.looped.bed -a $output/thisIN.bed3 > $output/thisIN.count.bed3
  # count the number of looped enhancers that in the largest possible encompassing IN region
  intersectBed -c  -F 1 -b $output/pOtherEnds.looped.bed -a $output/theseINs.merged.bed3 > $output/theseINs.count.bed3
  looped_total=$(wc -l $output/pOtherEnds.looped.bed|cut -d " " -f 1)
  echo $thisProm > $output/line
  # use line3 to store # of looped enhancers in the same INss
  echo $looped_total > $output/line2
  if [ -s $output/thisIN.count.bed3 ];then
    cat $output/thisIN.count.bed3 > $output/line3
    #paste $output/line $output/thisIN.count.bed3 $output/line2 >> $output/pergene.IN.txt
    looped_inIN=$(cut -f 4 $output/thisIN.count.bed3)
  else
    echo -e ".\t.\t.\t0" > $output/line3
    #paste $output/line $output/line3  $output/line2 >> $output/pergene.IN.txt
    looped_inIN=0
  fi
  # use line4 to store  # of looped enhancers in the same largest eompassing INs
  if [ -s $output/theseINs.count.bed3 ];then
    cat $output/theseINs.count.bed3 > $output/line4
    looped_inMergedIN=$(cut -f 4 $output/theseINs.count.bed3)
  else
    echo -e ".\t.\t.\t0" > $output/line4
    looped_inMergedIN=0
  fi


  paste $output/line $output/line3 $output/line4 $output/line2 >> $output/pergene.IN.txt
### SUMMARIZE PERGENE TRAITS
  if [ -s $output/thisIN.bed3 ]
  then
    export isIN="inIN"
  else
    export isIN="notIN"
  fi

  if [ -s $output/pOtherEnds.inIN.bed ]
  then
    export hasINenh="hasINenh"
  else
    export hasINenh="noINenh"
  fi

  if [ -s $output/pOtherEnds.looped.bed ]
  then
    export hasLoopEnh="hasLoopEnh"
  else
    export hasLoopEnh="noLoopEnh"
  fi

  echo -e "$thisProm\t$isIN,$hasINenh,$hasLoopEnh,$looped_total,$looped_inMergedIN,$looped_inIN,$pOtherEnds_enhancer,$pOtherEnds_StringentEnhancer,$pOtherEnds_promoter,$pOtherEnds_H3K4me3,$pOtherEnds_CTCF,$pOtherEnds_H3K27me3" >> $output/summary.txt
  #echo -e "$thisProm\t$isIN,$hasINenh,$hasLoopEnh,$looped_total,$looped_inMergedIN,$looped_inIN" >> $output/summary.txt


#### COUNT PROCESSED GENES SO FAR
  if ! (($counter % 100)); then
      echo "$counter promoters processed"
#      break
  fi
  counter=$((counter+1))

done < $output/prom.txt
#exit
echo -e "#######Find enhancers assignable multipe ways (loop and INs)#####\n"
#### FIND ENHS NOT ASSIGNABLE BY LOOPS OR NEIGHBORHOODS
cat $output/pergene.inIN.txt $output/pergene.looped.txt | awk -F\\t '{if($2 !="") print $2}' | sed s/\,/\\n/g | sed s/\_/\\t/g | awk -F\\t '{print $1 "\t" $2 "\t" $3}' > $output/AssignableRegions.bed
#intersectBed -v -a $enhancerSet -b $output/AssignableRegions.bed > $output/EnhancersNotAssignableByLoops.bed
intersectBed -v -a $output/theseEnhs.bed -b $output/AssignableRegions.bed > $output/EnhancersNotAssignableByLoops.bed

#### FIND ENHS OVERLAPPING GENES
intersectBed -loj -a $output/theseEnhs.bed -b $output/theseProms.bed | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/OverlappingProms.txt2
intersectBed -loj -a $output/theseEnhs.bed -b $wholegenes | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/OverlappingGenes.txt2
#intersectBed -loj -a $enhancerSet -b $promoters | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/OverlappingProms.txt2
#intersectBed -loj -a $enhancerSet -b $wholegenes | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/OverlappingGenes.txt2
#using new compound_lines.pl###
perl $ScriptDir/compound_lines.pl $output/OverlappingProms.txt2 0 1 > $output/OverlappingProms.txt
perl $ScriptDir/compound_lines.pl $output/OverlappingGenes.txt2 0 1 > $output/OverlappingGenes.txt

#### FIND NONLOOPASSIGNABLE OVERLAPPING GENES
intersectBed -loj -a $output/EnhancersNotAssignableByLoops.bed -b $promoters | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/NonLoopOverlappingProms.txt2
intersectBed -loj -a $output/EnhancersNotAssignableByLoops.bed -b $wholegenes | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/NonLoopOverlappingGenes.txt2
perl $ScriptDir/compound_lines.pl $output/NonLoopOverlappingProms.txt2 0 1 > $output/NonLoopOverlappingProms.txt
perl $ScriptDir/compound_lines.pl $output/NonLoopOverlappingGenes.txt2 0 1 > $output/NonLoopOverlappingGenes.txt


#### FIND ENHS BY PROXIMITY TO GENES
closestBed -t all -a $enhancerSet -b $promoters | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/Closest.txt2
perl $ScriptDir/compound_lines.pl $output/Closest.txt2 0 1 > $output/Closest.txt
#awk -F\\t '{if($2 !="") print $2}' $output/Closest.txt| sed s/\,/\\n/g | sed s/\_/\\t/g | awk -F\\t '{print $1 "\t" $2 "\t" $3}'|sort -k1,1 -k2,2n > $output/Closest.bed

#### FIND ENHS NOT ASSIGNABLE BY LOOPS OR OVERLAP
cat $output/pergene.inIN.txt $output/pergene.looped.txt $output/NonLoopOverlappingProms.txt $output/NonLoopOverlappingGenes.txt | awk -F\\t '{if($2 !="") print $2}' | sed s/\,/\\n/g | sed s/\_/\\t/g | awk -F\\t '{print $1 "\t" $2 "\t" $3}' > $output/AssignableRegions.bed

intersectBed -v -a $enhancerSet -b $output/AssignableRegions.bed | sort -k1,1 -k2,2n >| $output/EnhancersNotAssignableByLoopsOverlap.bed
closestBed -t all -a $output/EnhancersNotAssignableByLoopsOverlap.bed -b $promoters | awk -F\\t '{if($5>0) print $7 "\t" $1 "_" $2 "_" $3}' >| $output/RemainderClosest.txt2
perl $ScriptDir/compound_lines.pl $output/RemainderClosest.txt2 0 1 > $output/RemainderClosest.txt

#### CLEANUP
cd $output
rm -rf  line line2 line3 line4  pOtherEnds.inIN.bed pOtherEnds.inIN.bed2 pOtherEnds.looped.all.bed pOtherEnds.looped.all.bed2 pOtherEnds.looped.bed pOtherEnds.looped.enhancer.bed pOtherEnds.looped.enhancer.bed2 pOtherEnds.looped.noEnhancer.bed pOtherEnds.looped.noEnhancer.bed2 pOtherEnds.nonProm.inIN.bed AssignableRegions.bed *.txt2 *bed3 thisProm.bed
#rm -rf  line line2 line3 line4 pOtherEnds.CTCF.bed pOtherEnds.CTCF.bed2 pOtherEnds.enhancer.bed pOtherEnds.enhancer.bed2 pOtherEnds.H3K27me3.bed pOtherEnds.H3K27me3.bed2 pOtherEnds.H3K4me3.bed pOtherEnds.H3K4me3.bed2 pOtherEnds.inIN.bed pOtherEnds.inIN.bed2 pOtherEnds.looped.all.bed pOtherEnds.looped.all.bed2 pOtherEnds.looped.bed pOtherEnds.looped.enhancer.bed pOtherEnds.looped.enhancer.bed2 pOtherEnds.looped.noEnhancer.bed pOtherEnds.looped.noEnhancer.bed2 pOtherEnds.nonProm.inIN.bed AssignableRegions.bed *.txt2 *bed3 thisProm.bed
