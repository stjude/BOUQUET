Introduction
========================
A pipeline which takes input of
1)promoter and gene annotation
2)Enahncer set
3)Isoluted Neibourghoodloop set
4)Enhancer signal measurement(ChIP-Seq bam)

to

1) assign (3D) enhancers to its target genes according to chromatin loops (Loops) and isulated neighborhoods(INs)
2) establish enhancer/promoter interaction networks
3) define gene regulatory landscape through community detection method, directly connected neighborhooh or daisy-chained neighborhood.
4) call 3DSE associated genes by ranking per-gene total enhancer signal

Other utility includes:
1) assign linear enhancers to its target genes and further call Superenahcer associated genes by ranking the culmulative total enhancer signal per-gene (modified ROSE)
2) summarise the distribution of promoters and enhancers across different type of EP assignment evidence (eg. by loop, INs, overlap, etc...) *currently disabled*
3) compare SE genes difined by 3D and linear method *currently disabled*
4) annotate loop anchor with potential biological function (e.g., enhancer, promoter, isulator, etc...) *currently disabled*

author: Jie Lu and Brian Abraham
Abraham Lab
Department of Computational Biology
St. Jude Children’s Research Hospital


Change logs
=============================
v28 10/24/2023

1. Add TAD-limited linear proximity-based enhancer rescue
2. Fix double couting problem of promoter-overlaping enhancers
3. Add edge categery as a edge feature


v19 12/21/2022
---------------------
1. There is a fundament changes about how we define gene regulatory elements. Instead of E and P, now we use “bedtools merge” to get a set of non-overlapping CREs regions with E.label and P.label to indicate the identify of each region, in order to take into account the fact that some region might have both enhancer and promoter function.

2. For any analysis involving enhancer definition, except community-level enhancer signal loading calculation, enhancers are defined as CREs with exclusive enhancer identify  (with E.label===1 and P.label==0).

3. Promoters are defined as CREs as long as they are marked with P.label (P.label==1)
4. In the case of community-level enhancer signal loading calculation, CREs with both E.label and P.label are also considered as enhancers. And the community-level CREs signal loading calculation will no longer double-count region labeled as both enhancers and promoters.

5. For promoters shared by multiple isoforms/genes, only one representative "isoform|gene" id is selected to attach to it across analysis. A separate file “XXX.map” is generated to keep the mapping between each unused genes and its representative gene.

6. The parameters remains the same. Extra scripts are added to handle the process of re-defining enhancer/promoter set, recording collapsed genes with same promoters, switching definition of enhancer set according to different analysis etc.

=============================
v2 01/04/2022
---------------------
More detailed per-gene data are reported in the result table
1. community.enhancer
   The location of all enhancers included in each gene regulatory network
2. community.gene
   The gene ids included in each gene regulatory network, in the format of “Gene symbol|RefseqID”
3. community.promoter
   The location of all promoter regions (+/- 2KB of TSS) included in each gene regulatory network


Usage
========================
usage : bash 3DSE.wraper.sh -o OutDir [-d ScriptDir] [-p Promoters] [-I INs] [-l Loops] [-e Enhancers] [-s EnSignal] [-c EnControl] [-g HighlightGenes] [-h]
Use option -h for more information
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


Example usage based on a test dataset
=========================

cd /rgs01/project_space/abrahgrp/Baker_DIPG_CRC/common/3DSE/out

#Submiting pipeline to LSF

for S in 1;do ID=mES.mm10.newMED1.Ends1000.weightedLP.v28.LoopAndINRajWoGmFilterLarge.EnhOnly.FullEnhancer.PEAnchor.LowerPeakCutOff.woGm.K27acTwoThirds.RevisedLoop.TAD2012Plus2017Plus2023.$S;bsub -P "3DSE" -q compbio -J $ID -R "rusage[mem=70GB]" -N -oo $ID.o -eo $ID.e "bash /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script/3DSE_wrapper.v28.sh -o /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/out_jie/$ID -d /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script -p /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.txt -w /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/new_mm10.clean.wholegenes.bed.TPM1.filtered.woGm.K27acTwoThirds.txt -I /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mES_INs.mm10.Raj.woGm.sorted.FilterLarge.bed -l /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/HiChIP_mES_C3_UT_H3K27ac.Qvalue.0.01.3.0.4.H3K27ac_input_Peaks_1e-9.mm10.new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.combined.5000.20000000.bedpe -e /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K27ac_input_Peaks_1e-9.mm10.bed -s /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_MED1.sorted.bam -c /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_Input.sorted.bam -t /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/TADs2012Plus2017Plus2023.bed -g 'Irf2bpl|NM_145836.2,Sox2|NM_011443.4,Msc|NM_001360810.1,D7Ertd143e|NR_028425.1' -n -S $S -O";done


********* Please check the following inputs for the details of required format*********
********** Please note, all path should be provided as absolute path!****************
################Input papameters#####################
OutDir: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/out_jie/mES.mm10.newMED1.Ends1000.weightedLP.v28.LoopAndINRajWoGmFilterLarge.EnhOnly.FullEnhancer.PEAnchor.LowerPeakCutOff.woGm.K27acTwoThirds.RevisedLoop.TAD2012Plus2017Plus2023.1
ID: mES.mm10.newMED1.Ends1000.weightedLP.v28.LoopAndINRajWoGmFilterLarge.EnhOnly.FullEnhancer.PEAnchor.LowerPeakCutOff.woGm.K27acTwoThirds.RevisedLoop.TAD2012Plus2017Plus2023.1
ScriptDir: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script
Promoters: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.txt
WholeGenes: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/new_mm10.clean.wholegenes.bed.TPM1.filtered.woGm.K27acTwoThirds.txt
INs: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mES_INs.mm10.Raj.woGm.sorted.FilterLarge.bed
Loops: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/HiChIP_mES_C3_UT_H3K27ac.Qvalue.0.01.3.0.4.H3K27ac_input_Peaks_1e-9.mm10.new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.combined.5000.20000000.bedpe
Enhancers: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K27ac_input_Peaks_1e-9.mm10.bed
Enhancer Signal: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_MED1.sorted.bam
Enhancer Signal Control: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/Kagey_Input.sorted.bam
Highlight Genes: Irf2bpl|NM_145836.2,Sox2|NM_011443.4,Msc|NM_001360810.1,D7Ertd143e|NR_028425.1
flag for network figures: 1
flag for functional analysis of promoter-interacting regions: 1
flag for Loop Only edgess: 0
flag for allowing P to P assignment: 0
H3K27ac Peaks: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K27ac_input_Peaks_1e-5.mm10.bed
H3K4me3 Peaks: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K4me3_input_Peaks_1e-5.mm10.bed
CTCF Peaks: /research_jude/rgs01_jude/groups/abrahgrp/projects/3D_GENOME_CONSORTIUM/abrahgrp/Raj/INS_mESC/YoungLabChipseqmm10/YounglabData_CTCF_mm10_Seaseq_run-p9_kd-auto_peaks.bed
H3K27me3 Peaks: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/H3K27me3_input_Peaks_1e-3.mm10.bed
expression file: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/mm10.expressed.gene.txt
Seed value for R random process: 1
TADs: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/Jie/Generalizable_Code_v2/ref/TADs2012Plus2017Plus2023.bed






