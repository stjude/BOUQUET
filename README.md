<p align="center">

  <h1 align="center">
    [BOUQUET]
  </h1>

  <p align="center">
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template" target="_blank">
     <img alt="Status"
          src="https://img.shields.io/badge/status-active-success.svg" />
   </a>
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjudecloud/bioinformatics-tool-template"  />
   </a>
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/pulls"  target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjudecloud/bioinformatics-tool-template"  />
  </p>


  <p align="center">
   [Building Optimized Units of Quantified Enhancer Topologies] 
   <br />
   <a href="#"><strong>Explore the docs »</strong></a>
   <br />
   <a href="#"><strong>Read the paper »</strong></a>
   <br />
   <br />
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    | 
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
   <br />
    ⭐ Consider starring the repo! ⭐
   <br />
  </p>
</p>

---
## Authors
SOFTWARE AUTHORS: Jie Lv, Brian J Abraham
### Contact
jie.lu@stjude.org
 <br />
brian.abraham@stjude.org  


## Quick Start
Installing BOUQUET:

git clone https://github.com/stjude/BOUQUET.git; cd BOUQUET

## Description:
BOUQUET is an integrative, graph-theory-based approach that uses
multiple aspects of genome topology to probe communities of CREs, their bound apparatus, and their
target genes. It is a community-detection algorithm that uses concepts from graph theory to convert undirected networks
of CREs and their connections (“nodes” and “edges”, respectively) into discrete communities of CREs
associated with each expressed gene. 

Briefly, it is based on a modified label-propagation strategy wherein highly interconnected sets of nodes are unbiasedly
assembled into communities. Any additional CREs connected directly to a given promoter either by high-confidence HiChIP loops 
or shared INs are then incorporated into that gene’s community. Finally, any leftover CREs are added to a gene's community by their concurrence within the same TAD.


## Usage  

usage : bash 3DSE.wraper.sh -o OutDir [-d ScriptDir] [-p Promoters] [-I INs] [-l Loops] [-e Enhancers] [-s EnSignal] [-c EnControl] [-g HighlightGenes] [-h]
<br />
Use option -h for more information
<br />
    Options:
    <br />
        -o  OutDir              user specified directory for all output of the pipeline, recommend to reflect the unique combination of input files
        <br />
        -d  ScriptDir           directory for all 3DSE scripts
         <br />
        -p  Promoters           list of gene promoters, in the format of Chr, Start, End, TranscriptID|GeneID, default to be 4Kb region around TSS
        <br />
        -w  WholeGenes          list of whole gene body, in the format of Chr, Start, End, TranscriptID|GeneID, default to be whole gene body
        <br />
        -I  INs                 Isolated Neighborhoods, in the format of Chr, Start, End.
        <br />
        -l  Loops               Direct loops connecting two regulatory regions, in the format of bedpe. Default to be the ones defined from H3K27ac HiChIP.
        <br />
        -e  Enhancers           Set of enhancers in the format of Chr, Start, End. Default to be the ones defined from H3K27ac ChIP-Seq.
        <br />
        -s  EnSignal            ChiP-seq signal to quantify Enhancer signal intensity in the number of reads. The default is to be the Bam file of MED1 ChIP sample.
        <br />
        -c  EnControl           Control for ChiP-seq signal to quantify Enhancer signal intensity in the number of reads. Default to be Bam file of MED1 Control(Input) sample.
        <br />
        -g  HighlightGenes      a list of genes for which comprehensive GRN figures and tracks are to be generated. Expecting format 'GeneSymbol1|RefSeqID1,GeneSymbol2|RefSeqID2...'. Please note a single quote is needed because of the special character "|".
        <br />
        -n  Network             a switch indicating whether to draw network figures (7) for highlighted genes. could be slow for genes with complex landscapes. Default to turn off.
        <br />
        -O  OtherEnd            a switch indicating whether to perform functional characterization for promoter-interacting regions.
        <br />
        -f  LO                  a switch indicating whether to use both loops and INs-supported regulatory edges (default, 0) or only loop-supported regulatory edges (1).
        <br />
        -F  PromIN              a switch indicating whether to only assign E to P (default, 0) or assign both E and P to P (1) according to INs (default, 0).
        <br />
        -v  H3K27acPeaks        H3K27ac peaks to define enhancers
        <br />
        -x  H3K4me3Peaks        H3K4me3 peaks to define promoters
        <br />
        -y  CTCFPeaks           CTCF peaks to define insulators
        <br />
        -z  H3K27me3Peaks       H3K27me3 peaks to define PcG sites
        <br />
        -E  ExpFile             gene symbol-level expression file calculated based on RNA-Seq
        <br />
        -S  Seed                Seed number of R random process. same seed value would  generate reproducible gene communities based on the LP method.
        <br />
        -t  TAD                 Topological Associated Domain used for rescue EO (Enhancer Orphan) after community-based Enhancer assignment
        <br />


