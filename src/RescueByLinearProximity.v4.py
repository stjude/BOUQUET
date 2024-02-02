from __future__ import print_function
import pybedtools
from collections import defaultdict
import time
import sys
import os.path as op
import argparse
import math
import commands
###V.1 updates
### 1. report distance to closet promter for each community enhancer and average distance to promoter for all community enhancers (showcase the unique ablility of our framework to discover distal enhancer-promoter connection)
### 2. reprot distance of each community enhancer to its focal promoter/gene and average distance of all community enhancers, giving a picture of how a gene's contributing regulatory elements distribute in terms of linear genomes.

###V.2 updates
###add a extra column to update table, community.CREs(community.E + community.P)

###V.3 updates
###add a extra column to update table, exp (expression level calculated based on PolyA RNA-seq)

usage = (
"""
    Rescue unassigned enhancers and combine them with community assigned results. Only enhancers rescused by method 1 are incorporated for downstrem analysis.

    Progresive linear proximity-based enhancer assignment:
    1. Within 10kb around promoter(4kb region)(since the resolution of HiChIP is 5kb)
    2. Within 10kb around gene body
    3. assign the rest based on linear promximity to promoters

"""
)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
def load_exp(file):
    """
    store exp level for each gene
    """
    start_time=time.time()
    if file:
        f=open(file)
    gene2exp=defaultdict(float)
    while True:
        line=f.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        gene,exp=cols
        exp=float(exp)
        gene2exp[gene]=exp
    f.close()
    eprint("time cost for load gene expression: ",time.time()-start_time, "seconds")
    return gene2exp

def locate_enhancer(file):
    """
    get clostest gene and the distance of each enhancer to its clostest gene
    """
    start_time=time.time()
    if file:
        f=open(file)
    enh2prom=defaultdict(list)
    enhDist=defaultdict(int)
    while True:
        line=f.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        enh,enh_id,gene,dist,signal=cols
        dist, signal=map (int, (dist,signal))
        enh2prom[enh].append(gene)
        enhDist[enh]=dist
    f.close()
    eprint("time cost for locate_enhancer: ",time.time()-start_time, "seconds")
    return enh2prom,enhDist


def map_enhancer(file,cutoff):
    """
    map enhancers based on linear proximity and distance cutoff
    """
    eprint("enhancer linear assignment file: ", file)
    start_time=time.time()
    if file:
        f=open(file)
    gene2enh=defaultdict(list)
    gene2enhLoad=defaultdict(int)
    while True:
        line=f.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        enh,enh_id,gene,dist,signal=cols
        dist, signal=map (int, (dist,signal))
        if dist <= cutoff:
            gene2enh[gene].append(enh)
            gene2enhLoad[gene]+=signal
    f.close()
    eprint("time cost for map_enhancer: ",time.time()-start_time, "seconds")
    return gene2enh,gene2enhLoad
def cal_dist(e,p):
    """
    calculate the absolute distance between given enhancer and promoter.

    """
    e_chr,e_s,e_e=e.split('_')
    p_chr,p_s,p_e=p.split('_')

    if e_chr != p_chr:
        eprint("ERROR!! Promoter and enhancer are not in the same chromosome: ",p,e)
        sys.exit(1)
    e_s,e_e,p_s,p_e=map (int, (e_s,e_e,p_s,p_e))
    if e_e<p_s or p_e<e_s:
    # cases for seperate:
    #1. enhancer to the left of promoter
    #2. promoter to the left of enhancer
        dist=abs((e_s+0.5*(e_e-e_s))-(p_s+0.5*(p_e-p_s)))
    else:
    #case for overlap
        dist=0
    return dist


def update_network(file,outfile,gene2enh_prom,gene2enhLoad_prom,gene2enh_genebody,gene2enhLoad_genebody,gene2enh_rest,gene2enhLoad_rest,enh2prom,enhDist,gene2exp):
    """
    given the community output file, update the file with extra columns related to linear rescued enhancers, enhancer distance to clostest genes, and enhancer distance to focal gene
    """
    eprint("network file: ", file)
    start_time=time.time()
    if file:
        f=open(file)
        of=open(outfile,"w")
    while True:
        line=f.readline()
        if line.startswith('id\t'):
            of.write(line.strip()+"\t"+"\t".join(("prom.proximity.enhancer","prom.proximity.enhancer.signal","genebody.proximity.enhancer","genebody.proximity.enhancer.signal","rest.proximity.enhancers","rest.proximity.enhancers.signal","total.E.signal","essential.total.E.signal","enhancer.closestGene","enhancer.distToClosestGene","enhancer.distToFocalGene","enhancer.AverageDistToClosestGene","enhancer.AverageDistToFocalGene","community.order","exp"))+"\n")
            continue
        if not line:
            break
        cols=line.strip().split("\t")
        gene=cols[0]
        gene_id,tx_id=gene.split("|")
        exp=gene2exp[gene_id]
        med1_comm=int(cols[30])
        community_order=int(cols[24])+int(cols[25])
        enhs=cols[32]
        genes=cols[33]
        promLoc=cols[34]
        promFocal=promLoc.split(",")[genes.split(",").index(gene)]
        #print (gene,promFocal)
        distToClosest=[]
        closestGene=[]
        distToFocal=[]
        if enhs != "":
            for e in enhs.split(","):
                distToClosest.append(enhDist[e])
                closestGene.append(";".join(enh2prom[e]))
                distToFocal.append(cal_dist(e,promFocal))
            average_distToClosest=sum(distToClosest)/float(len(distToClosest))
            average_distToFocal=sum(distToFocal)/float(len(distToFocal))
        else:
             distToClosest=['NA']
             closestGene=["NA"]
             distToFocal=["NA"]
             average_distToClosest="NA"
             average_distToFocal="NA"

        if gene not in gene2enh_prom:
            gene2enh_prom[gene]=["NA"]
            gene2enhLoad_prom[gene]=0
        if gene not in gene2enh_genebody:
            gene2enh_genebody[gene]=["NA"]
            gene2enhLoad_genebody[gene]=0
        if gene not in gene2enh_rest:
            gene2enh_rest[gene]=["NA"]
            gene2enhLoad_rest[gene]=0
        total_signal=med1_comm+gene2enhLoad_prom[gene]+gene2enhLoad_genebody[gene]+gene2enhLoad_rest[gene]
        essential_total_signal=med1_comm+gene2enhLoad_prom[gene]
        of.write(line.strip()+"\t"+"\t".join(map (str, (",".join(gene2enh_prom[gene]),gene2enhLoad_prom[gene],",".join(gene2enh_genebody[gene]),gene2enhLoad_genebody[gene],",".join(gene2enh_rest[gene]),gene2enhLoad_rest[gene],total_signal,essential_total_signal,",".join(closestGene),",".join(map(str,distToClosest)),",".join(map(str,distToFocal)),average_distToClosest,average_distToFocal,community_order,exp)))+"\n")
    f.close()
    of.close()
    #for gene in gene2exp:
    #    print (gene, gene2exp[gene])
    eprint("time cost for update_network: ",time.time()-start_time, "seconds")

def main():
    """
    Rescue community missed enhancers and combine them with community assigned results. Only enhancers rescused by promoter neighbour region are incorporated for downstrem analysis.
    """
    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument(
        "-n",
        "--network",
        dest="network",
        help="output table of community assignment. extra column of 1)location of linear proximity rescued enhancers, 2) med1 signal of all assigned enhancers (community + promoter neighbour rescued) "
    )
    ap.add_argument(
        "-a",
        "--promP",
        dest="promProx",
        help="enhancers missed by community assigned by promoter proximity"
    )
    ap.add_argument(
        "-b",
        "--GeneBodyP",
        dest="geneBodyProx",
        help="enhancers missed by promoter proximity assignment assigned by genebody proximity"
    )
    ap.add_argument(
        "-c",
        "--rest",
        dest="restPromProx",
        help="enhancers missed by genebody proximity assignment (the rest) assigned by promoter proximity"
    )
    ap.add_argument(
        "-e",
        "--enhs",
        dest="enhs",
        help="txt file containing enhancers mapped to closest genes with distance information, output from bedtool clostest function"
    )
    ap.add_argument(
        "-E",
        "--exp",
        dest="exp",
        help="txt file containing gene-level expression value calculated from PolyA RNA-seq"
    )

    args = ap.parse_args()

    if (not args.network) or (not args.promProx) or ( not args.geneBodyProx) or (not args.restPromProx) or (not args.enhs) or (not args.exp):
        ap.print_help()
        sys.exit(1)
    gene2enh_prom,gene2enhLoad_prom=map_enhancer(args.promProx,5000)
    gene2enh_genebody,gene2enhLoad_genebody=map_enhancer(args.geneBodyProx, 5000)
    gene2enh_rest,gene2enhLoad_rest=map_enhancer(args.restPromProx, float('inf'))
    enh2prom,enhDist=locate_enhancer(args.enhs)
    gene2exp=load_exp(args.exp)
    outfile=args.network+".update"
    update_network(args.network,outfile,gene2enh_prom,gene2enhLoad_prom,gene2enh_genebody,gene2enhLoad_genebody,gene2enh_rest,gene2enhLoad_rest,enh2prom,enhDist,gene2exp)


    #for gene in gene2enh_rest:
    #    print(gene,gene2enhLoad_rest[gene],gene2enh_rest[gene])

if __name__ == "__main__":
    main()

