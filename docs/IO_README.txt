######
script
######
3DSE_wrapper.v28.sh

######
input files
######
ID: mES.mm10.newMED1.Ends1000.weightedLP.v28.LoopAndINRajWoGmFilterLarge.EnhOnly.FullEnhancer.PEAnchor.LowerPeakCutOff.woGm.K27acTwoThirds.RevisedLoop.TAD2012Merge.CTCFUpdate.1
ScriptDir: /research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/script
Promoters: new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.txt
WholeGenes: new_mm10.clean.wholegenes.bed.TPM1.filtered.woGm.K27acTwoThirds.txt
INs: mES_INs.mm10.Raj.woGm.sorted.FilterLarge.bed (need Raj's input)
Loops: HiChIP_mES_C3_UT_H3K27ac.Qvalue.0.01.3.0.4.H3K27ac_input_Peaks_1e-9.mm10.new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.combined.5000.20000000.bedpe  (GSM2645434,GSM2645435,GSM2645436)
TAD: MergedDixon2012_mm10.bed
Enhancers: H3K27ac_input_Peaks_1e-9.mm10.bed (GSM1526287)
MED1 Signal: Kagey_MED1.sorted.bam (GSM560347,GSM560348; need Brian's double check )
MED1 Signal Control: Kagey_Input.sorted.bam (GSM560357; need Brian's double check )
H3K27ac Peaks: H3K27ac_input_Peaks_1e-5.mm10.bed (GSM1526287)
H3K4me3 Peaks: H3K4me3_input_Peaks_1e-5.mm10.bed (GSM1526288)
CTCF Peaks: YounglabData_CTCF_mm10_Seaseq_run-p9_kd-auto_peaks.bed (SRR299029)
H3K27me3 Peaks: H3K27me3_input_Peaks_1e-3.mm10.bed (GSM1526291)
H3K27ac Bam: SRR1613251.bam
H3K4me3 Bam: SRR1613252.bam
CTCF Bam: SRR299029.bam
H3K27me3 Bam: SRR1613255.bam
expression file: mm10.expressed.gene.txt (GSM1246869)
Seed value for R random process: 1
flag for network figures: 1
flag for functional analysis of promoter-interacting regions: 1
flag for Loop Only edgess: 0
flag for allowing P to P assignment: 0

######
Output files
######
OutDir:/research_jude/rgs01_jude/groups/abrahgrp/projects/Baker_DIPG_CRC/common/3DSE/out_jie/mES.mm10.newMED1.Ends1000.weightedLP.v28.LoopAndINRajWoGmFilterLarge.EnhOnly.FullEnhancer.PEAnchor.LowerPeakCutOff.woGm.K27acTwoThirds.RevisedLoop.TAD2012Merge.CTCFUpdate.1

"all-in-one" community table:
network.plusIN.promoter.11050.seed.1.PE_nodes.Kagey_MED1.sorted.bam.txt.update
network.plusIN.promoter.11050.seed.1.PE_nodes.Kagey_MED1.sorted.bam.txt.update.uniqueCommunity

######
utility files
######

Parsimony INs assigned: INs.parsimony.bed
enhancers without overlapping promoters: H3K27ac_input_Peaks_1e-9.mm10.bed.new
unique promoters: new_mm10.clean.4kbproms.bed.TPM1.filtered.woGm.K27acTwoThirds.txt.new
mm10_genome: the folder contains the mm10 RefSeq ref genes used in this project, check the README for details
gene_list: BOUQUET and ROSE defiend SE genes. For BOUQUET, four gene lists are generated: 1)E only without TAD rescue; 2) E+P without TAD rescue; 3) E only with TAD rescue; 4) E+P with TAD rescue. Please note all genes with their pomoters collapsed in network building are uncolloped now.
