##new in 4.0 version
#1. output network.P and network. E in the table
#2. generate modified network figures for publication.

## new in 5.0 version
#print chromsome region correspoinding to the entire network of focal gene

## new in 6.0 
# set seed as a parameter instead of hard coded to evaluate the stability of gene community inferred by weighted LP community detection method. 

## new in 7.0
#1. Use less greedy community detection method: instead of assigning focal gene and its direct neighbors the same initial label, use default LP method (each node has its unique label) with direct neighbors back-filling.
#2. report unique community ID for each gene: componentID.communityMembershipID

## new in 7.3
#1. new ways to visulize community
#2. distinguish edges as loop only, IN only and loop+IN supported 

## new in 7.4
#1. add ability to calculating CRE signal loading by substracting promoter signal

##new in 7.6

# now community.E is non-promoter Enhancers, but  the enhancer list contains all enhancers, whether or not overlap with promoters
############################
## weightedLP with backfilling,rescure any direct neighbours that are not included in the final community 
library(igraph)
library(dplyr)
# caution: must set the seed (set.seed(1)) before any random draw and sample to get reproducible results.
args = commandArgs(trailingOnly=TRUE)
#wd="/Users/jlu/Project/HiChIP/3DSE/mm10/LP_test_20230123"
#n_file="PE_nodes.Kagey_MED1.sorted.bam.txt"
#e_file="PE_edges.txt"
#genes='H4c6|NM_175655.3,H3c6|NM_178205.3,H2ac8|NM_178187.4,H2bc4|NM_001290380.1,Irf2bpl|NM_145836.2,Msc|NM_001360810.1,H2ac10|NM_175661.3,Suz12|NM_001163018.1,Sox2|NM_011443.4,Sulf2|NM_001252578.1,D7Ertd143e|NR_028425.1,Pou5f1|NM_013633.3,Nanog|NM_001289831.1,Klf4|NM_010637.3,H2-M5|NM_001115075.2,Gsk3a|NM_001031667.1,Cenpt|NM_177150.2,Mmachc|NM_025962.3'
#i='Msc|NM_001360810.1'
#i='H2ac10|NM_175661.3'
#i='H2ac8|NM_178187.4'
#i='H2bc4|NM_001290380.1'
#args=c(wd,n_file,e_file,genes)
# test if there is at least 3 arguments: if not, return an error
if (length(args)<3) {
    stop("At least three argument must be supplied (working directory, nodes file, and edge file).\n", call.=FALSE)
} else {
    # print args
    working_dir<-args[1]
    node_file<-args[2]
    edge_file<-args[3]
    highlightGenes<-c()
    seed<-1
    cat(paste0("Working directory: ",working_dir,"\n","Node file: ", node_file,"\n","Edge file: ", edge_file,"\n"))
    if (length(args)>=4) {
        cat(paste0("Highlighted genes: ", args[4],"\n"))
        highlightGenes<-strsplit(args[4],',',fixed=TRUE)[[1]]
    }
    if (length(args)==5) {
        cat(paste0("seed value: ", args[5],"\n"))
        seed<-strtoi(args[5])
    }

}
    
setwd(working_dir)


nodes <- read.csv(node_file, header=T, as.is=T, sep = "\t")
links <- read.csv(edge_file, header=T, as.is=T,sep = "\t")

cat(paste0("Number of nodes: ", nrow(nodes)),"\n")
cat(paste0("Number of unique nodes: ", length(unique(nodes$id))),"\n")
cat(paste0("Number of edges: ", nrow(links)),"\n")
cat(paste0("Number of edges of the same type (either Loop or IN supported) with unique start and end point: ",nrow(unique(links[,c("from", "to","type")])),"\n"))
cat(paste0("Number of total edges with unique start and end point: ",nrow(unique(links[,c("from", "to")])),"\n"))

#only keep edges supported by loop
#links<-links[links$type==0,]

# when wanting to merge edges of same origin (either loop or IN supported) 
#links <- aggregate(links[,4], links[,-4], sum)
#links <- links[order(links$from, links$to),]
#colnames(links)[4] <- "weight"
#rownames(links) <- NULL

#build the graph
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
#plot(net, edge.arrow.size=.4, vertex.label=NA)

cat(paste0("Number of nodes before simplify: ", vcount(net)),"\n")
cat(paste0("Number of edges before simplify: ", ecount(net)),"\n")
# remove self connected edges, and merge edges with same start and end point, if edges are supported by both Loop and IN, mark the edge with stronger evidence, i.e., Loop (type=0 for loop and type=1 for IN) 
#net <- simplify(net, remove.multiple = T, remove.loops = T,edge.attr.comb=list(weight="max",type="min")) 
#instead of mark the edge with a single stronger evidence (loop if both supported by loop and IN), use one of the 3 types (0 for loop, 1 for IN, and 01 for Loop+IN)

net <- simplify(net, remove.multiple = T, remove.loops = T,edge.attr.comb=list(weight="max",type=toString)) 
# mearge multiple assignment e.g. from "11010" to "01" 
E(net)$type<-sapply(E(net)$type,function(x) paste(sort(unique(strsplit(x,', ')[[1]])),collapse = ''))
## 2 is the edge supported by both loops and INs
E(net)$type[E(net)$type=='01']<-2
cat(paste0("Number of nodes after simplify: ", vcount(net)),"\n")
cat(paste0("Number of edges after simplify: ", ecount(net)),"\n")

#calculate proportion of edges from INs 
tab<-table(E(net)$type)
cat(paste0("proportion of edges supported by Loop only : ",tab[1]," ",sum(tab)," ",tab[1]/(sum(tab))),"\n")
cat(paste0("proportion of edges supported by IN only : ",tab[2]," ",sum(tab)," ",tab[2]/(sum(tab))),"\n")
cat(paste0("proportion of edges supported by Loop+IN: ",tab[3]," ",sum(tab)," ",tab[3]/(sum(tab))),"\n")

#Loop through all Promoters

V(net)$name_full=paste(V(net)$name1,V(net)$name2,sep='|')


comp=components(net)
eigen=eigen_centrality(net)
V(net)$comp=comp$membership
V(net)$comp_size=comp$csize[V(net)$comp]
V(net)$degree=degree(net)
V(net)$size=log(V(net)$degree)
V(net)$eigen=eigen$vector

#Sort promoter by degree
order2=order(V(net)$degree[V(net)$Plabel==1])

#vlist<-head(V(net)$name_full[V(net)$group==1][rev(order2)],n=30)
#vertex_attr(net, index = V(net)[name_full%in%vlist])

#Sort promoter by degree
promoters=V(net)$name_full[V(net)$Plabel==1][rev(order2)]

N<-length(promoters)

#Output edgelist and edge type ( 0="Loop only"; 1="IN only"; 2="Both")
el<-as_edgelist(net)
edges<-data.frame(node1=el[,1],node2=el[,2],type=E(net)$type)
edge_file<-paste0("network.plusIN.promoter.",N,".",node_file,".edges")
write.table(edges, edge_file, row.names=FALSE, col.names=TRUE, sep="\t", quote=F)

# all output in "network" sub-directory
outdir<-file.path(working_dir,'network.LP.default')
dir.create(outdir)
setwd(outdir)

output<-paste0("network.plusIN.promoter.",N,".seed.",seed,".",node_file)
cat(paste0("Output file: ", output,"\n"))

if (! file_test("-f",output)) {
    cat(paste0("Performing gnome-wide network inference and generating output: ", output,"\n"))

    #Dynamically growing dataframe in an efficient way

    DF <- data.frame(id=rep("", N), 
                  membership=rep(NA, N),
                  degree=rep(NA, N),
                  network.order=rep(NA,N),
                  network.P=rep(NA,N),
                  network.E=rep(NA,N),
                  network.size=rep("",N),
                  network.density=rep(NA,N),
                  neighborhood.E.order1=rep(NA,N),
                  neighborhood.P.order1=rep(NA,N),
                  neighborhood.size.order1=rep(NA,N),
                  neighborhood.density.order1=rep(NA,N),
                  neighborhood.interaction.order1=rep(NA,N),
                  neighborhood.interaction.density.order1=rep(NA,N),
                  neighborhood.E.signal.order1=rep(NA,N),
                  neighborhood.signal.order1=rep(NA,N),
                  neighborhood.E.order2=rep(NA,N),
                  neighborhood.P.order2=rep(NA,N),
                  neighborhood.size.order2=rep(NA,N),
                  neighborhood.density.order2=rep(NA,N),
                  neighborhood.interaction.order2=rep(NA,N),
                  neighborhood.interaction.density.order2=rep(NA,N),
                  neighborhood.E.signal.order2=rep(NA,N),
                  neighborhood.signal.order2=rep(NA,N),
                  community.E=rep(NA,N),
                  community.P=rep(NA,N),
                  community.size=rep(NA,N),
                  community.density=rep(NA,N),
                  community.interaction=rep(NA,N),
                  community.interaction.density=rep(NA,N),
                  community.E.signal=rep(NA,N),
                  community.signal=rep(NA,N),
                  community.enhancer=rep("",N),
                  community.gene=rep("",N),
                  community.promoter=rep("",N),
                  stringsAsFactors=FALSE) 
    # create progress bar
    pb <- txtProgressBar(min = 0, max = N, style = 3)

    for (i in 1:N){
    # caution: must set the seed (set.seed(1)) before any random draw and sample to get reproducible results.
    set.seed(seed)
        # focal promoter
        p=promoters[i]
        # network of p
        landscape=subcomponent(net,V(net)[V(net)$name_full==p])
        graph_i<-induced_subgraph(net,landscape)
  
        # get distance of all other nodes to node i
        order<-distances(graph_i, v=V(graph_i)[name_full==p], weights=NA)
  
        #total enhancer signal of order 1 neighbors for node i
        neighborhood.signal.order1<-sum(V(graph_i)[order==1]$Signal)
  
        #total enhancer signal of order 1 neighbors for node i, only count neighbors that are enhancers
        neighborhood.E.signal.order1<-sum(V(graph_i)[order==1][Elabel==1]$SignalSub)
  
        #total enhancer signal of order 1 and 2 neighbors for node i
        neighborhood.signal.order2<-sum(V(graph_i)[order==1|order==2]$Signal)
  
        #total enhancer signal of order 1 and 2 neighbors for node i, only count neighbors that are enhancers
        neighborhood.E.signal.order2<-sum(V(graph_i)[order==1|order==2][Elabel==1]$SignalSub)
  
  
        # get order1 and order2 neighbors
        neighbors.nodes<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==p], mindist=1)[[1]]
        neighbors.nodes2<-ego(graph_i,order=2, nodes=V(graph_i)[V(graph_i)$name_full==p], mindist=2)[[1]]
  
        #get full neighborhood edges
        neighbors.graph<-make_ego_graph(graph_i,order=1, nodes=V(graph_i)[name_full==p], mindist=0)[[1]]
        neighbors.graph2<-make_ego_graph(graph_i,order=2, nodes=V(graph_i)[name_full==p], mindist=0)[[1]]
  
        # get neighborhood interaction edges only
        neighbors.graph3<-make_ego_graph(graph_i,order=1, nodes=V(graph_i)[name_full==p], mindist=1)[[1]]
        neighbors.graph4<-make_ego_graph(graph_i,order=2, nodes=V(graph_i)[name_full==p], mindist=2)[[1]]
  
        #get edges of shortest paths to order 1 neighbors
        path.s1<-lapply(all_shortest_paths(graph_i,V(graph_i)[name_full==p],to=neighbors.nodes,weights=NA)$res,as_ids)
        path.s1<-sapply(X=path.s1,FUN=get.edge.ids,graph=graph_i)
  
        #get edges of shortest paths to both order 1 and 2 neighbors
        #this is a workaround to get edge ids of all shortest paths from p to order2(N) neighbors
        path.s2<-lapply(all_shortest_paths(graph_i,V(graph_i)[name_full==p],to=union(neighbors.nodes,neighbors.nodes2),weights=NA)$res,function(x) unlist(combn(x,2,simplify=F)))
        path.s2<-unlist(sapply(X=path.s2,FUN=get.edge.ids,graph=graph_i))
        path.s2<-unique(path.s2[path.s2!=0])
   
        #get all the summary statistics
        membership<-V(net)$comp[V(net)$name_full==p]
        degree<-V(net)$degree[V(net)$name_full==p]
        network.order<-length(landscape)
        #Output edgelist and edge type ( 0="Loop only"; 1="IN only"; 2="Both")
        network.size<-paste(ecount(graph_i),length(E(graph_i)[E(graph_i)$type=="0"]),length(E(graph_i)[E(graph_i)$type=="1"]),length(E(graph_i)[E(graph_i)$type=="2"]),sep=',')
        network.density<-graph.density(graph_i)
        network.P<-sum(V(graph_i)$Plabel==1)
        network.E<-sum(V(graph_i)$Elabel==1)
        network.E2<-sum(V(graph_i)$Elabel==1 & V(graph_i)$Plabel==0)
        network.E<-paste(network.E,network.E2,sep=',')
  
        #E and P number of neighbors
        P<-sum(neighbors.nodes$Plabel==1)
        E<-sum(neighbors.nodes$Elabel==1)
  
        P2<-sum(neighbors.nodes2$Plabel==1)
        E2<-sum(neighbors.nodes2$Elabel==1)
  
        neighborhood.size.order1<-ecount(neighbors.graph)
        neighborhood.size.order2<-ecount(neighbors.graph2)
  
        neighborhood.interaction.order1<-ecount(neighbors.graph3)
        neighborhood.interaction.order2<-ecount(neighbors.graph4)
  
        neighborhood.density.order1<-graph.density(neighbors.graph)
        neighborhood.density.order2<-graph.density(neighbors.graph2)
  
        neighborhood.interaction.density.order1<-graph.density(neighbors.graph3)
        neighborhood.interaction.density.order2<-graph.density(neighbors.graph4)
  
        #Community detection based on edge betweenness 
        #ceb <- cluster_edge_betweenness(graph_i,weights=NULL)
  
        #Community detection based on greedy optimization of modularity
        #cfg <- cluster_fast_greedy(graph_i)
  
        #Community detection based on propagating labels
        #clp <- cluster_label_prop(graph_i,weights=NULL)
  
        #m<-membership(clp)[V(graph_i)[name_full==p]]
        #community.nodes<-communities(clp)[[m]]
  
        #the subgraph of entire community 
        #community<-induced_subgraph(graph_i,community.nodes) 
        
        # the subgraph after removing i
        #community.landscape<- induced_subgraph(graph_i,community.nodes[!community.nodes==V(graph_i)[name_full==p]$name])
        
        #community.P<-sum(V(community)$group==1)
        #community.E<-sum(V(community)$group==0)
        #community.size<-ecount(community)
        #community.density<-graph.density(community)
        #community.interaction<-ecount(community.landscape)
        #community.interaction.density<-graph.density(community.landscape)
        #community.E.signal<-sum(V(community.landscape)[group.label=='E']$Signal)
        #community.signal<-sum(V(community.landscape)$Signal)
        
        #######Community detection based on weighted propagating labels
        # find all direct neighbors of a focal gene
        #neighbors.nodes.all<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==p], mindist=0)[[1]]
        # Initialize all direct neighbors with the same label (the same label as focal gene), with the rest nodes of the network keeping their own unique label (a way to take into acount the fact that direct neighbors are the most strong candidates for a given gene's regulators)   
        #initial.status<-ifelse(as.numeric(V(graph_i)) %in% as.numeric(V(graph_i)[neighbors.nodes.all]),as.numeric(V(graph_i)[name_full==p])-1,as.numeric(V(graph_i))-1)
        #fix the label of focal promoter and its direct neighbors, meaning forcelly including all direct neighours in the final community
        #fixed<-initial.status==as.numeric(V(graph_i)[name_full==p])-1 
        #clp2 <- cluster_label_prop(graph_i,weights = NULL,initial=initial.status,fixed=fixed)
        #backing filling
        #clp2 <- cluster_label_prop(graph_i,weights = NULL,initial=initial.status)
        clp <- cluster_label_prop(graph_i,weights = NA)

        #m2<-membership(clp2)[V(graph_i)[name_full==p]]
        m<-membership(clp)[V(graph_i)[name_full==p]]
        # instead of only report componentID, resport componentID.communityMembershipID, so genes within the same community will have the same unique ID. 
        membership<-paste0(membership,".",m)
        #community.nodes2<-communities(clp2)[[m2]]
        community.nodes<-communities(clp)[[m]]
        ##back filling, rescure any direct neighbours that are not included in the final community
        #community.nodes2<- unique(c(community.nodes2,names(neighbors.nodes)))
        community.nodes2<- unique(c(community.nodes,names(neighbors.nodes)))

        community.nodes2.enhancer<-paste(as_ids(V(graph_i)[community.nodes2][Elabel==1]),collapse=",")
        community.nodes2.full_name<-paste(V(graph_i)[community.nodes2][Plabel==1]$name_full,collapse=",")
        community.nodes2.promoter<-paste(as_ids(V(graph_i)[community.nodes2][Plabel==1]),collapse=",")

        #the subgraph of entire community
        community2<-induced_subgraph(graph_i,community.nodes2)
        # the subgraph after removing i
        community.landscape2<- induced_subgraph(graph_i,community.nodes2[!community.nodes2==V(graph_i)[name_full==p]$name])

        #network.E<-sum(V(graph_i)$Elabel==1)
        #network.E2<-sum(V(graph_i)$Elabel==1 & V(graph_i)$Plabel==0)
        #network.E<-paste(network.E,network.E2,sep=',')
        community.P2<-sum(V(community2)$Plabel==1)
        community.E2<-sum(V(community2)$Elabel==1)
        community.E3<-sum(V(community2)$Elabel==1 & V(community2)$Plabel==0)
        community.E2<-paste(community.E2,community.E3,sep=',')

        community.size2<-paste(ecount(community2),length(E(community2)[E(community2)$type=="0"]),length(E(community2)[E(community2)$type=="1"]),length(E(community2)[E(community2)$type=="2"]),sep=',')
        #network.size<-paste(ecount(graph_i),length(E(graph_i)[E(graph_i)$type=="0"]),length(E(graph_i)[E(graph_i)$type=="1"]),length(E(graph_i)[E(graph_i)$type=="2"]),sep=',')
        community.density2<-graph.density(community2)
        community.interaction2<-ecount(community.landscape2)
        community.interaction.density2<-graph.density(community.landscape2)
        #community.E.signal2<-sum(V(community.landscape2)[group.label=='E']$Signal)
        #community.signal2<-sum(V(community.landscape2)$Signal)  
        community.E.signal2<-sum(V(community2)[Elabel==1]$SignalSub)
        #when measuring enhancer loading using promoter+enhancer, the focal promoter itself should also be counted 
        community.signal2<-sum(V(community2)$Signal)  
        
        #print(c(p,community.E,community.P,community.size,community.density,community.interaction,community.interaction.density,community.E.signal,community.signal))
  
        DF[i,]<-list(p,membership,degree,network.order,network.P,network.E,network.size,network.density,
                    E,P,neighborhood.size.order1,neighborhood.density.order1,neighborhood.interaction.order1,neighborhood.interaction.density.order1,neighborhood.E.signal.order1,neighborhood.signal.order1,
                    E2,P2,neighborhood.size.order2,neighborhood.density.order2,neighborhood.interaction.order2,neighborhood.interaction.density.order2,neighborhood.E.signal.order2,neighborhood.signal.order2,
                    community.E2,community.P2,community.size2,community.density2,community.interaction2,community.interaction.density2,community.E.signal2,community.signal2,community.nodes2.enhancer,community.nodes2.full_name,community.nodes2.promoter)
        # update progress bar
        setTxtProgressBar(pb, i)
  
    } 

    close(pb)

    write.table(DF, output, row.names=FALSE, col.names=TRUE, sep="\t", quote=F)
}else {
    cat(paste0("The following output file detected, skipping gnome-wide network inference: ", output,"\n"))
}

#example gene list
#highlightGenes<-c("Irf2bpl|NM_145836","Ahsa1|NM_146036","Tmed8|NM_001033475","Cipc|NM_173735","Esrrb|NM_001159500","Esrrb|NM_011934","Suz12|NM_199196","Hist1h2bh|NM_178197","Hist1h4h|NM_153173")
# 'Irf2bpl|NM_145836,Suz12|NM_199196,Hist1h4h|NM_153173'
#Irf2bpl|NM_145836.2,Msc|NM_001360810.1
#i='Msc|NM_001360810.1'
#i='H2ac10|NM_175661.3'
#i='H2ac8|NM_178187.4'
#i='H2bc4|NM_001290380.1'

if (length(highlightGenes)>0) {

for (i in highlightGenes){
    set.seed(seed)
    ##transform from geneID|txID to txID|geneID
    #geneID=unlist(strsplit(i,'|',fixed=TRUE))[1]
    #txID=unlist(strsplit(i,'|',fixed=TRUE))[2]
    #i<-paste0(txID,"|",geneID)
    
    landscape=subcomponent(net,V(net)[V(net)$name_full==i])
    graph_i<-induced_subgraph(net,landscape)
  
  
    #Community detection based on edge betweenness 
    ceb <- cluster_edge_betweenness(graph_i,weights=NULL)

    #Community detection based on greedy optimization of modularity
    cfg <- cluster_fast_greedy(graph_i,weights=NULL)

    #Community detection based on propagating labels
    clp <- cluster_label_prop(graph_i,weights = NULL)
    # find all direct neighbors of a focal gene
    neighbors.nodes.all<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==i], mindist=0)[[1]]
    # Initialize all direct neighbors with the same label (the same label as focal gene), with the rest nodes of the network keeping their own unique label (a way to take into acount the fact that direct neighbors are the most strong candidates for a given gene's regulators)   
    initial.status<-ifelse(as.numeric(V(graph_i)) %in% as.numeric(V(graph_i)[neighbors.nodes.all]),as.numeric(V(graph_i)[name_full==i])-1,as.numeric(V(graph_i))-1)
    #fix the label of focal promoter and its direct neighbors, meaning forcelly including all direct neighours in the final community
    #fixed<-initial.status==as.numeric(V(graph_i)[name_full==i])-1 
    set.seed(seed)
    #clp2 <- cluster_label_prop(graph_i,weights = NULL,initial=initial.status,fixed=fixed)
    clp2 <- cluster_label_prop(graph_i,weights = NULL,initial=initial.status)
    
    neighbors.nodes<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==i], mindist=1)[[1]]
    ######## Default LP #########
    
    m<-membership(clp)[V(graph_i)[name_full==i]]
    community.nodes<-communities(clp)[[m]]
  
    #the subgraph of entire community 
    community<-induced_subgraph(graph_i,community.nodes) 
    # the subgraph after removing i
    community.landscape<- induced_subgraph(graph_i,community.nodes[!community.nodes==V(graph_i)[name_full==i]$name])
  
    community.P<-sum(V(community)$Plabel==1)
    community.E<-sum(V(community)$Elabel==1)
    community.size<-ecount(community)
    community.density<-graph.density(community)
    community.interaction<-ecount(community.landscape)
    community.interaction.density<-graph.density(community.landscape)
    community.E.signal<-sum(V(community.landscape)[Elabel==1]$SignalSub)
    community.signal<-sum(V(community.landscape)$Signal)
    
    ######## Default LP with backfilling#########
    
    community.nodes3<- unique(c(community.nodes,names(neighbors.nodes)))
    community.rescued.nodes<-setdiff(names(neighbors.nodes),community.nodes)  
    community.nodes3.enhancer<-as_ids(V(graph_i)[community.nodes3][Elabel==1])%>% gsub('_',"\t",.)
    community.nodes3.full_name<-V(graph_i)[community.nodes3][Plabel==1]$name_full
    community.nodes3.promoter<-as_ids(V(graph_i)[community.nodes3][Plabel==1])%>% gsub('_',"\t",.)
    df.community.nodes.genes<-data.frame(community.nodes3.full_name,community.nodes3.promoter,stringsAsFactors=FALSE)
    if (length(V(graph_i)[community.rescued.nodes][Elabel==1])>0) {
        community.rescued.enhancer<-as_ids(V(graph_i)[community.rescued.nodes][Elabel==1])%>% gsub('_',"\t",.)
        #replace "|" by ".", since the "|" in the file name give proteinpaint error message. 
        file.community<-paste0(i,".seed.",seed,".community.rescued.enhancer.LPDefault.bed") %>% gsub('\\|','.',.)
        fileConn<-file(file.community)
        writeLines(community.rescued.enhancer%>% gsub('_',"\t",.),fileConn)
        close(fileConn)
    }
    if (length(V(graph_i)[community.rescued.nodes][Plabel==1])>0) {
        community.rescued.promoter<-as_ids(V(graph_i)[community.rescued.nodes][Plabel==1])%>% gsub('_',"\t",.)
        community.rescued.gene<-V(graph_i)[community.rescued.nodes][Plabel==1]$name_full
        df.community.rescued.genes<-data.frame(community.rescued.gene,community.rescued.promoter,stringsAsFactors=FALSE)
        write.table(df.community.rescued.genes, paste0(i,".seed.",seed,".community.rescued.genes.LPDefault.bed")%>% gsub('\\|','.',.), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
    }
    
    #output gene list of community
    write.table(df.community.nodes.genes, paste0(i,".seed.",seed,".community.genes.LPDefault.bed")%>% gsub('\\|','.',.), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
  
  
    #replace "|" by ".", since the "|" in the file name give proteinpaint error message. 
    file.community<-paste0(i,".seed.",seed,".community.enhancer.LPDefault.bed") %>% gsub('\\|','.',.)
    fileConn<-file(file.community)
  
    #output enhancer list of community
    writeLines(community.nodes3.enhancer%>% gsub('_',"\t",.),fileConn)
    close(fileConn)
  
    #output the linear chromosomal range of community
  
    community.chr<-strsplit(community.nodes3[1],"_")[[1]][1]
    community.start<-strsplit(community.nodes3[1],"_")[[1]][2]
    community.end<-strsplit(community.nodes3[length(community.nodes3)],"_")[[1]][3]
    fileConn2<-file(paste0(i,".seed.",seed,".community.range.LPDefault.bed")%>% gsub('\\|','.',.))
    writeLines(c(community.chr,community.start,community.end),sep='\t',fileConn2)
    close(fileConn2)
    
    #the subgraph of entire community 
    community3<-induced_subgraph(graph_i,community.nodes3) 
    # the subgraph after removing i
    community.landscape3<- induced_subgraph(graph_i,community.nodes3[!community.nodes3==V(graph_i)[name_full==i]$name])
  
    community.P3<-sum(V(community3)$Plabel==1)
    community.E3<-sum(V(community3)$Elabel==1)
    community.size3<-ecount(community3)
    community.density3<-graph.density(community3)
    community.interaction3<-ecount(community.landscape3)
    community.interaction.density3<-graph.density(community.landscape3)
    community.E.signal3<-sum(V(community3)[Elabel==1]$SignalSub)
    community.signal3<-sum(V(community3)$Signal)

    
    ######## Weighted LP with backfilling#########
  
    m2<-membership(clp2)[V(graph_i)[name_full==i]]
    community.nodes2<-communities(clp2)[[m2]]
    ##back filling, rescure any direct neighbours that are not included in the final community
    #neighbors.nodes<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==i], mindist=1)[[1]]
    community.nodes2<- unique(c(community.nodes2,names(neighbors.nodes)))
    community.rescued.nodes2<-setdiff(names(neighbors.nodes),community.nodes2)  
  
    community.nodes2.enhancer<-as_ids(V(graph_i)[community.nodes2][Elabel==1])%>% gsub('_',"\t",.)
    community.nodes2.full_name<-V(graph_i)[community.nodes2][Plabel==1]$name_full
    community.nodes2.promoter<-as_ids(V(graph_i)[community.nodes2][Plabel==1])%>% gsub('_',"\t",.)
  
    df.community.nodes.genes<-data.frame(community.nodes2.full_name,community.nodes2.promoter,stringsAsFactors=FALSE)
    if (length(V(graph_i)[community.rescued.nodes2][Elabel==1])>0) {
        community.rescued2.enhancer<-as_ids(V(graph_i)[community.rescued.nodes2][Elabel==1])%>% gsub('_',"\t",.)
        #replace "|" by ".", since the "|" in the file name give proteinpaint error message. 
        file.community<-paste0(i,".seed.",seed,".community.rescued.enhancer.bed") %>% gsub('\\|','.',.)
        fileConn<-file(file.community)
        writeLines(community.rescued2.enhancer%>% gsub('_',"\t",.),fileConn)
        close(fileConn)
    }
    if (length(V(graph_i)[community.rescued.nodes2][Plabel==1])>0) {
        community.rescued2.promoter<-as_ids(V(graph_i)[community.rescued.nodes2][Plabel==1])%>% gsub('_',"\t",.)
        community.rescued2.gene<-V(graph_i)[community.rescued.nodes2][Plabel==1]$name_full
        df.community.rescued2.genes<-data.frame(community.rescued2.gene,community.rescued2.promoter,stringsAsFactors=FALSE)
        write.table(df.community.rescued2.genes, paste0(i,".seed.",seed,".community.rescued.genes.bed")%>% gsub('\\|','.',.), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
    }
    
    #output gene list of community
    write.table(df.community.nodes.genes, paste0(i,".seed.",seed,".community.genes.bed")%>% gsub('\\|','.',.), row.names=FALSE, col.names=FALSE, sep="\t", quote=F)
  
  
    #replace "|" by ".", since the "|" in the file name give proteinpaint error message. 
    file.community<-paste0(i,".seed.",seed,".community.enhancer",".bed") %>% gsub('\\|','.',.)
    fileConn<-file(file.community)
  
    #output enhancer list of community
    writeLines(community.nodes2.enhancer%>% gsub('_',"\t",.),fileConn)
    close(fileConn)
  
    #output the linear chromosomal range of community
  
    community.chr<-strsplit(community.nodes2[1],"_")[[1]][1]
    community.start<-strsplit(community.nodes2[1],"_")[[1]][2]
    community.end<-strsplit(community.nodes2[length(community.nodes2)],"_")[[1]][3]
    fileConn2<-file(paste0(i,".seed.",seed,".community.range.bed")%>% gsub('\\|','.',.))
    writeLines(c(community.chr,community.start,community.end),sep='\t',fileConn2)
    close(fileConn2)
    
    #output the linear chromosomal range of network
    network.nodes<-attr(landscape,which="names")
    network.chr<-strsplit(network.nodes[order(sapply(strsplit(network.nodes, "_"), function(x) x[2]))][1],"_")[[1]][1]
    network.start<-strsplit(network.nodes[order(sapply(strsplit(network.nodes, "_"), function(x) x[2]))][1],"_")[[1]][2]
    network.end<-strsplit(network.nodes[order(sapply(strsplit(network.nodes, "_"), function(x) x[3]))][length(network.nodes)],"_")[[1]][3]


    fileConn3<-file(paste0(i,".seed.",seed,".network.range.bed")%>% gsub('\\|','.',.))
    writeLines(c(network.chr,network.start,network.end),sep='\t',fileConn3)
    close(fileConn3)

    #the subgraph of entire community 
    community2<-induced_subgraph(graph_i,community.nodes2) 
    # the subgraph after removing i
    community.landscape2<- induced_subgraph(graph_i,community.nodes2[!community.nodes2==V(graph_i)[name_full==i]$name])
  
    community.P2<-sum(V(community2)$Plabel==1)
    community.E2<-sum(V(community2)$Elabel==1)
    community.size2<-ecount(community2)
    community.density2<-graph.density(community2)
    community.interaction2<-ecount(community.landscape2)
    community.interaction.density2<-graph.density(community.landscape2)
    community.E.signal2<-sum(V(community2)[Elabel==1]$SignalSub)
    community.signal2<-sum(V(community2)$Signal)
  
    cat(paste0( "GeneSymbol|RefSeqID: ",i,"\n",
              "################## Community  based on default Label Propagation #############","\n",
              "Number of Community Enhancer: ",community.E,"\n",
              "Number of Community Promoter: ",community.P,"\n",
              "Number of Community Edges: ", community.size,"\n",
              "Community Density: ", community.density, "\n",
              "Signal loading for Community Enhancer: ", community.E.signal,"\n",
              "Signal loading for Community Enhancer plus Promoter: ",community.signal, "\n\n"))
    cat(paste0( "################## Community based on default Label Propagation with backfilling#############","\n",
              "Number of Community Enhancer: ",community.E3,"\n",
              "Number of Community Promoter: ",community.P3,"\n",
              "Number of Community Edges: ", community.size3,"\n",
              "Community Density: ", community.density3, "\n",
              "Signal loading for Community Enhancer: ", community.E.signal3,"\n",
              "Signal loading for Community Enhancer plus Promoter: ",community.signal3, "\n\n"))
    cat(paste0( "################## Community based on weighted Label Propagation with backfilling#############","\n",
              "Number of Community Enhancer: ",community.E2,"\n",
              "Number of Community Promoter: ",community.P2,"\n",
              "Number of Community Edges: ", community.size2,"\n",
              "Community Density: ", community.density2, "\n",
              "Signal loading for Community Enhancer: ", community.E.signal2,"\n",
              "Signal loading for Community Enhancer plus Promoter: ",community.signal2, "\n\n"))

    # get distance of all other nodes to node i
    order<-distances(graph_i, v=V(graph_i)[name_full==i], weights=NA)
  
    #total enhancer signal of order 1 neighbors for node i
    neighborhood.signal.order1<-sum(V(graph_i)[order==1]$Signal)
  
    #total enhancer signal of order 1 neighbors for node i, only count neighbors that are enhancers
    neighborhood.E.signal.order1<-sum(V(graph_i)[order==1][Elabel==1]$SignalSub)
  
    #total enhancer signal of order 1 and 2 neighbors for node i
    neighborhood.signal.order2<-sum(V(graph_i)[order==1|order==2]$Signal)
  
    #total enhancer signal of order 1 and 2 neighbors for node i, only count neighbors that are enhancers
    neighborhood.E.signal.order2<-sum(V(graph_i)[order==1|order==2][Elabel==1]$SignalSub)
  
    # get order1 and order2 neighbors
   # neighbors.nodes<-ego(graph_i,order=1, nodes=V(graph_i)[V(graph_i)$name_full==i], mindist=1)[[1]]
    neighbors.nodes2<-ego(graph_i,order=2, nodes=V(graph_i)[V(graph_i)$name_full==i], mindist=2)[[1]]
    
    #get full neighborhood edges
    neighbors.graph<-make_ego_graph(graph_i,order=1, nodes=V(graph_i)[name_full==i], mindist=0)[[1]]
    neighbors.graph2<-make_ego_graph(graph_i,order=2, nodes=V(graph_i)[name_full==i], mindist=0)[[1]]
  
    # get neighborhood interaction edges only
    neighbors.graph3<-make_ego_graph(graph_i,order=1, nodes=V(graph_i)[name_full==i], mindist=1)[[1]]
    neighbors.graph4<-make_ego_graph(graph_i,order=2, nodes=V(graph_i)[name_full==i], mindist=2)[[1]]
  
  
    #get edges of shortest paths to order 1 neighbors
    path.s1<-lapply(all_shortest_paths(graph_i,V(graph_i)[name_full==i],to=neighbors.nodes,weights=NA)$res,as_ids)
    path.s1<-sapply(X=path.s1,FUN=get.edge.ids,graph=graph_i)
  
    #get edges of shortest paths to both order 1 and 2 neighbors
    path.s2<-lapply(all_shortest_paths(graph_i,V(graph_i)[name_full==i],to=union(neighbors.nodes,neighbors.nodes2),weights=NA)$res,function(x) unlist(combn(x,2,simplify=F)))
    path.s2<-unlist(sapply(X=path.s2,FUN=get.edge.ids,graph=graph_i))
    path.s2<-unique(path.s2[path.s2!=0])
    
    #set up edge color
    ecol <- rep("gray80", ecount(graph_i))
    ecol[path.s2] <- "blue"
    ecol[path.s1] <- "tomato"
    E(graph_i)$edge.color<-ecol
    edge.evidence.color<-c('pink',"SkyBlue2","red")
    E(graph_i)$edge.evidence.color<-edge.evidence.color[as.numeric(E(graph_i)$type)+1]
  
    #find neighbors of focal promoter
    P<-sum(neighbors.nodes$Plabel==1)
    E<-sum(neighbors.nodes$Elabel==1)
  
    P2<-sum(neighbors.nodes2$Plabel==1)
    E2<-sum(neighbors.nodes2$Elabel==1)
  
    membership<-V(net)$comp[V(net)$name_full==i]
    #componentID.communityID
    membership<-paste0(membership,".",m)
    comp.size<-V(net)$comp_size[V(net)$name_full==i]
    #summary statistics
    cat(paste0( "################## Gene Regulatory Landscape defined by first and second order neighborhood #############","\n",
              "Unique network ID: ",membership,"\n",
              "Size of network: ",comp.size,"\n",
              "Degree of gene Promoter: ", V(net)$degree[V(net)$name_full==i],"\n",
              "Number of first-order neiborhood Enhancer (directly connected): ",E, "\n",
              "Number of first-order neiborhood Promoter (directly connected): ",P, "\n",
              "Number of second-order neiborhood Enhancer (daisy-chained): ",E2, "\n",
              "Number of second-order neiborhood Promoter (daisy-chained): ",P2, "\n\n\n"))
    #set up vertex color
  
    #colrs <- c("lightsteelblue2", "#555555")
    #colrs <- c("lightsteelblue2", "black")
    #vcol<-colrs[V(graph_i)$Plabel+1]
    #vcol[V(graph_i)$name_full==i] <- "tomato"
    #V(graph_i)$color<-vcol
  
    # determine layout 
    set.seed(seed)
    l <- layout_with_fr(graph_i) 
   
#set up vertex color
  
    colrs <- c("SkyBlue2","orange")
    shapes<-c("square","circle")
    sizes<-c(3,4)
    
    
    vshape<-shapes[V(graph_i)$Plabel+1]
    vcol<-colrs[V(graph_i)$Plabel+1]
    vsize<-sizes[V(graph_i)$Plabel+1]
    #vcol[V(graph_i)$name_full==i] <- "tomato"
    #vcol[V(graph_i)$name %in% community.rescued.nodes] <- "yellow"
    
    V(graph_i)$color<-vcol
    V(graph_i)$shape<-vshape
    V(graph_i)$size<-vsize
    #vcol[V(graph_i)$name_full==i] <- "tomato"
   # V(graph_i)$color[V(graph_i)$name_full==i] <-'tomato'
    V(graph_i)$size[V(graph_i)$name_full==i] <-8
    
    E(graph_i)$width<-1
    E(graph_i)$width[path.s1]<-1
    
   # V(graph_i)$community<-clp$membership
    #rain <- rainbow(length(unique(V(graph_i)$community)), alpha=.5)
    
    #E(graph_i)$color<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community == V(graph_i)[x[2]]$community, rain[V(graph_i)[x[1]]$community], '#D3D3D3'))
    
    #plot(graph_i, vertex.size=4, vertex.label=NA, edge.color=E(graph_i)$color,layout=l)
    #set up vertex label display
  
    #V(graph_i)$label<-ifelse(V(graph_i)$Plabel==1, V(graph_i)$name1, NA)
    V(graph_i)$label<-ifelse(V(graph_i)$Plabel==1, V(graph_i)$name1, V(graph_i)$name)
    # separately draw direct and indirect edges of the component
  
    graph_i.rest1<-graph_i-neighbors.graph  
    graph_i.rest2<-graph_i-neighbors.graph2 
    graph_i.neighbor1<-graph_i-graph_i.rest1
    graph_i.neighbor2<-graph_i-graph_i.rest2
  
    #only highlight shortest paths instead of entire neighborhood network
    graph_i.p.rest1<-graph_i-E(graph_i)[path.s1]
    graph_i.p.rest2<-graph_i-E(graph_i)[path.s2]
    graph_i.p.neighbor1<-graph_i-graph_i.p.rest1
    graph_i.p.neighbor2<-graph_i-graph_i.p.rest2
    graph_i.rest3<-graph_i-neighbors.graph3
    graph_i.rest4<-graph_i-neighbors.graph4
    graph_i.neighbor3<-graph_i-graph_i.rest3
    graph_i.neighbor4<-graph_i-graph_i.rest4
  
  
############## Draw Fighures ####################### 

######################### Whole network #####################
#change the seperator between gene symbol and RefSeq ID from "|" to  "." in the output file.
  file1=paste0(i,".seed.",seed,".network.pdf") %>% gsub('\\|','.',.)
  pdf(file=file1)
  plot(graph_i,
       layout=l,
       edge.color=E(graph_i)$edge.color,
       edge.width=E(graph_i)$width,
       vertex.label.dist=1,
       vertex.label.cex=0.7,
       vertex.shapes=V(graph_i)$shape,
       vertex.label.color="black")
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("Network Order:", length(landscape)),
    paste("Network Promoters:", sum(V(graph_i)$Plabel==1)),
    paste("Network Enhancers:", sum(V(graph_i)$Elabel==1)),
    paste("Network Size:", ecount(graph_i)),
    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
    paste("Network Membership:", membership),
    paste("N1 Enhancers:",E ),
    paste("N1 Promoters:",P ),
    paste("N2 Enhancers:",E2 ),
    paste("N2 Promoters:",P2 ),
    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
    paste("N1 Signal:",neighborhood.signal.order1),
    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
    paste("N2 Signal:",neighborhood.signal.order2)
  ),
    cex=0.75
  )
 
  
  dev.off()
  
######################### Whole network with edge colored by evidence (orange:loop, blue:IN, red:Loop+IN) #####################
  file8=paste0(i,".seed.",seed,".network.edgeColoredByEvidence.pdf") %>% gsub('\\|','.',.)
  pdf(file=file8)
  plot(graph_i,
       layout=l,
       edge.color=E(graph_i)$edge.evidence.color,
       vertex.label.dist=0.6,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.shapes=V(graph_i)$shape)
  legend("topleft", legend=c(
    paste("Loop only:", length(E(graph_i)[E(graph_i)$type=="0"])),
    paste("IN only:", length(E(graph_i)[E(graph_i)$type=="1"])),
    paste("Loop and IN:", length(E(graph_i)[E(graph_i)$type=="2"]))
    #paste("Gene:",i),
    #paste("Network Order:", length(landscape)),
    #paste("Network Promoters:", sum(V(graph_i)$Plabel==1)),
    #paste("Network Enhancers:", sum(V(graph_i)$Elabel==1)),
    #paste("Network Size:", ecount(graph_i)),
    #paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i])
  ),
  cex=0.75
  )
  
  dev.off()
######################### Neighborhood Edges Comparison #####################
  
  file2=paste0(i,".seed.",seed,".neighborhood.edges.pdf") %>% gsub('\\|','.',.)
  pdf(file=file2,width=18,height=12)
  par(mfrow=c(1,3))
  plot(graph_i.p.rest2,
       layout=l,
       vertex.label.dist=1,
       edge.color=E(graph_i.p.rest2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Edge: rest")
  plot(graph_i.p.neighbor2,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.p.neighbor2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Edge: order2")
  plot(graph_i.p.neighbor1,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.p.neighbor1)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Edge: direct")
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("Network Order:", length(landscape)),
    paste("Network Size:", ecount(graph_i)),
    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
    paste("N1 Enhancers:",E ),
    paste("N1 Promoters:",P ),
    paste("N2 Enhancers:",E2 ),
    paste("N2 Promoters:",P2 ),
    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
    paste("N1 Signal:",neighborhood.signal.order1),
    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
    paste("N2 Signal:",neighborhood.signal.order2)
  ))
  dev.off()
  
######################### Neighborhood Density Comparison #####################
  file3=paste0(i,".seed.",seed,".neighborhood.DensityComparison.pdf")  %>% gsub('\\|','.',.)
  pdf(file=file3,width=18,height=12)
  par(mfrow=c(1,2))
  
  plot(graph_i.neighbor3,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.neighbor3)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Neighborhood Interaction: Order1")
  plot(graph_i.neighbor4,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.neighbor3)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Neighborhood Interaction: Order2")
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("N1 Interaction Density:",graph.density(neighbors.graph3)),
    paste("N2 Interaction Density:",graph.density(neighbors.graph4)),
    paste("N1 Interaction:",ecount(graph_i.neighbor3)),
    paste("N2 Interaction:",ecount(graph_i.neighbor4))
  ))
  dev.off()
  
######################### Neighborhood Subnetwork Comparison #####################
  file4=paste0(i,".seed.",seed,".neighborhood.subnetwork.pdf") %>% gsub('\\|','.',.)
  pdf(file=file4,width=18,height=12)
  par(mfrow=c(1,3))
  plot(graph_i.rest2,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.rest2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Other")
  plot(graph_i.neighbor2,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.neighbor2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Neighborhood: order2")
  plot(graph_i.neighbor1,
       layout=l,
       vertex.label.dist=0.6,
       edge.color=E(graph_i.neighbor1)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA,
       main="Neighborhood: order1")
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("Network Order:", length(landscape)),
    paste("Network Size:", ecount(graph_i)),
    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
    paste("N1 Enhancers:",E ),
    paste("N1 Promoters:",P ),
    paste("N2 Enhancers:",E2 ),
    paste("N2 Promoters:",P2 ),
    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
    paste("N1 Signal:",neighborhood.signal.order1),
    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
    paste("N2 Signal:",neighborhood.signal.order2)
  ))
  dev.off()
  
  # perform community detection

######################### Community Detection Methods Comparison #####################
  
  file5=paste0(i,".seed.",seed,".community.allMethods.pdf") %>% gsub('\\|','.',.)
  pdf(file=file5,width=18,height=12)
  par(mfrow=c(2,2))
  plot(ceb, graph_i,layout=l,vertex.label=NA,main="Community detection based on edge betweenness")
  
  plot(cfg, graph_i,layout=l,vertex.label=NA,main="Community detection based on greedy optimization of modularity")
  
  plot(clp, graph_i,layout=l, vertex.label=NA,main="Community detection based on propagating labels")
  
  plot(graph_i.p.neighbor2,
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.p.neighbor2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Direct edges")
  
  dev.off()
  ######################### Loop VS. IN Supported Edges Comparison #####################
  
  graph_i.Loop <- graph_i - E(graph_i)[E(graph_i)$type=="1" | E(graph_i)$type=="2" ]
  graph_i.IN <- graph_i - E(graph_i)[E(graph_i)$type=="0"| E(graph_i)$type=="2"]
  graph_i.loopAndIN<-graph_i - E(graph_i)[E(graph_i)$type=="0" | E(graph_i)$type=="1" ]
  
  file7=paste0(i,".seed.",seed,".LoopVsIN.edges.pdf")  %>% gsub('\\|','.',.)
  
  pdf(file=file7,width=18,height=12)
  
  par(mfrow=c(1,3))
  
  plot(graph_i.Loop,
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.Loop)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Edge: Loop supported")
  
  plot(graph_i.IN, 
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.IN)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Edge: IN supported")
  legend("topleft", legend=c(
    paste("Edge by Loop:",ecount(graph_i.Loop)),
    paste("Edge by IN:",ecount(graph_i.IN)),
    paste("Edge by LoopAndIN:",ecount(graph_i.loopAndIN))
  ))
  
  plot(graph_i.loopAndIN, 
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.loopAndIN)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Edge: Loop and IN supported")
  dev.off() 
  
  ##reset edge color to only color direct edges
  
  ecol <- rep("gray80", ecount(graph_i))
  #ecol[path.s2] <- "blue"
  ecol[path.s1] <- "tomato"
  E(graph_i)$edge.color<-ecol
  
  ######################### Whole network without label  #####################
  #A version of whole network without color for second degree nodes and promoter label
  file8=paste0(i,".seed.",seed,".network.noLabel.pdf") %>% gsub('\\|','.',.)
  pdf(file=file8)
  plot(graph_i,
       layout=l,
       edge.color=E(graph_i)$edge.color,
       vertex.label.dist=0.6,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       vertex.label=NA)
  
  dev.off()
  
######################### Community Based on Weighted Label Propagation Method #####################

file6=paste0(i,".seed.",seed,".weighted.label.propagation.community.pdf") %>% gsub('\\|','.',.)
  
  pdf(file=file6,width=18,height=12)
  par(mfrow=c(2,3))
  #plot(graph_i.p.neighbor2,
  #     layout=l,
  #     vertex.label.dist=0.6,
  #     vertex.label=NA,
  #     edge.color=E(graph_i.p.neighbor2)$edge.color,
  #     vertex.label.cex=0.7,
  #     vertex.label.color="black",
  #     main="Direct and Daisy-Chain Neighbors")
  
  
  plot(clp, graph_i,layout=l, vertex.label=NA,main="Default LP Community")
#  legend("bottomleft", legend=c(
#    paste("Gene:",i),
#    paste("Network Order:", length(landscape)),
#    paste("Network Size:", ecount(graph_i)),
#    paste("Network Membership:", membership),
#    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
#    paste("N1 Enhancers:",E ),
#    paste("N1 Promoters:",P ),
#    paste("N2 Enhancers:",E2 ),
#    paste("N2 Promoters:",P2 ),
#    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
#    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
#    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
#    paste("N1 Signal:",neighborhood.signal.order1),
#    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
#    paste("N2 Signal:",neighborhood.signal.order2)
#  ))
  op <- par(cex = 0.5)
  legend("bottomleft", legend=c(
    
    paste("Community Enhancers:",community.E ),
    paste("Community Promoters:",community.P ),
    paste("Community Size:",community.size ),
    paste("Community Enhancer Signal:",community.E.signal),
    paste("Community Signal:",community.signal)
  ))
  legend("bottom", legend=c(
    
    paste("Community Enhancers:",community.E3 ),
    paste("Community Promoters:",community.P3 ),
    paste("Community Size:",community.size3 ),
    paste("Community Enhancer Signal:",community.E.signal3),
    paste("Community Signal:",community.signal3)
  ))
  legend("bottomright", legend=c(
    
    paste("Community Enhancers:",community.E2 ),
    paste("Community Promoters:",community.P2 ),
    paste("Community Size:",community.size2 ),
    paste("Community Enhancer Signal:",community.E.signal2),
    paste("Community Signal:",community.signal2)
  ))
  
  
  #vcol[V(graph_i)$name %in% community.rescued.nodes] <- "yellow"
  #ecol[path.s1] <- "tomato"
  #E(graph_i)$edge.color<-ecol
  #V(graph_i)$color<-vcol
 
  V(graph_i)$community<-clp$membership
  V(graph_i)$community.revise<-clp$membership
  V(graph_i)$community.revise[V(graph_i)$name %in% community.rescued.nodes]<- membership(clp)[V(graph_i)$name_full==i]
  
  
  rain <- rainbow(max(V(graph_i)$community), alpha=.5)
  rain.revise<- rainbow(max(V(graph_i)$community.revise), alpha=.5)
  
  #V(graph_i)$color<-rain[membership(clp)]
  V(graph_i)$color<-rain[V(graph_i)$community]
  V(graph_i)$color.revise<-rain.revise[V(graph_i)$community.revise]
  
  #V(graph_i)$color[V(graph_i)$name_full==i] <-'tomato'
  V(graph_i)$color[V(graph_i)$name %in% community.rescued.nodes] <-rain[membership(clp)[V(graph_i)$name_full==i]]
  E(graph_i)$color<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community == V(graph_i)[x[2]]$community, '#454545', '#D3D3D3'))
  E(graph_i)$color[path.s1] <- "tomato"
  
  #V(graph_i)$color[V(graph_i)$name %in% community.rescued.nodes] <-rain[membership(clp)[V(graph_i)$name_full==i]]
  E(graph_i)$color.revise<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community.revise == V(graph_i)[x[2]]$community.revise, '#454545', '#D3D3D3'))
  E(graph_i)$color.revise[path.s1] <- "tomato"
  
  graph_i2<-graph_from_data_frame(d<-igraph::as_data_frame(graph_i,what="edges")%>%arrange(factor(color,levels=c("#D3D3D3","#454545","tomato"))),vertices=igraph::as_data_frame(graph_i,what="vertices"),directed = FALSE)
  graph_i2.revise<-graph_from_data_frame(d<-igraph::as_data_frame(graph_i,what="edges")%>%arrange(factor(color.revise,levels=c("#D3D3D3","#454545","tomato"))),vertices=igraph::as_data_frame(graph_i,what="vertices"),directed = FALSE)
  
  plot(graph_i2, edge.width=E(graph_i2)$width,vertex.frame.width=0.1,vertex.size=V(graph_i2)$size, vertex.color=V(graph_i)$color, vertex.label=NA,vertex.shapes=V(graph_i2)$shape,edge.color=E(graph_i2)$color,layout=l)
  plot(graph_i2.revise, edge.width=E(graph_i2.revise)$width,vertex.frame.width=0.1,vertex.size=V(graph_i2.revise)$size, vertex.color=V(graph_i)$color.revise,vertex.label=NA,vertex.shapes=V(graph_i2.revise)$shape,edge.color=E(graph_i2.revise)$color.revise,layout=l)
  
  #plot(graph_i, edge.width=E(graph_i)$width,vertex.size=V(graph_i)$size,vertex.frame.width=0.1, vertex.label=NA,vertex.shapes=V(graph_i)$shape, edge.color=E(graph_i)$color,layout=l)
  #set up vertex color
  
    #colrs <- c("lightsteelblue2","orange")
    #vcol<-colrs[V(graph_i)$Plabel+1]
    #vcol[V(graph_i)$name_full==i] <- "tomato"
    #vcol[V(graph_i)$name %in% community.rescued.nodes2] <- "yellow"
    #V(graph_i)$color<-vcol
    
    V(graph_i)$community2<-clp2$membership
    V(graph_i)$community2.revise<-clp2$membership
    V(graph_i)$community2.revise[V(graph_i)$name %in% community.rescued.nodes] <- membership(clp2)[V(graph_i)$name_full==i]
    
    #rain <- rainbow(length(unique(V(graph_i)$community2)), alpha=.5)
    rain <- rainbow(max(V(graph_i)$community2), alpha=.5)
    rain.revise<- rainbow(max(V(graph_i)$community2.revise), alpha=.5)
    
    #V(graph_i)$color2<-rain[membership(clp2)]
    V(graph_i)$color2<-rain[V(graph_i)$community2]
    V(graph_i)$color2.revise<-rain.revise[V(graph_i)$community2.revise]
    #V(graph_i)$color[V(graph_i)$name_full==i] <-'tomato'
    
    V(graph_i)$color2[V(graph_i)$name %in% community.rescued.nodes] <- rain[membership(clp2)[V(graph_i)$name_full==i]]
    
    #E(graph_i)$color<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community == V(graph_i)[x[2]]$community, rain[V(graph_i)[x[1]]$community], '#D3D3D3'))
    E(graph_i)$color2<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community2 == V(graph_i)[x[2]]$community2, '#454545', '#D3D3D3'))
    E(graph_i)$color2[path.s1] <- "tomato"
    
    E(graph_i)$color2.revise<-apply(as.data.frame(get.edgelist(graph_i)), 1, function(x) ifelse(V(graph_i)[x[1]]$community2.revise == V(graph_i)[x[2]]$community2.revise, '#454545', '#D3D3D3'))
    E(graph_i)$color2.revise[path.s1] <- "tomato"
  
    plot(clp2, graph_i,layout=l, vertex.label=NA,main="Weighted LP Community")
    graph_i2<-graph_from_data_frame(d<-igraph::as_data_frame(graph_i,what="edges")%>%arrange(factor(color2,levels=c("#D3D3D3","#454545","tomato"))),vertices=igraph::as_data_frame(graph_i,what="vertices"),directed = FALSE)
    graph_i2.revise<-graph_from_data_frame(d<-igraph::as_data_frame(graph_i,what="edges")%>%arrange(factor(color2.revise,levels=c("#D3D3D3","#454545","tomato"))),vertices=igraph::as_data_frame(graph_i,what="vertices"),directed = FALSE)
    
    plot(graph_i2, edge.width=E(graph_i2)$width,vertex.frame.width=0.1,vertex.size=V(graph_i2)$size, vertex.color=V(graph_i)$color2,vertex.label=NA,vertex.shapes=V(graph_i2)$shape,edge.color=E(graph_i2)$color2,layout=l)
    plot(graph_i2.revise, edge.width=E(graph_i2.revise)$width,vertex.frame.width=0.1,vertex.size=V(graph_i2.revise)$size, vertex.color=V(graph_i)$color2.revise,vertex.label=NA,vertex.shapes=V(graph_i2.revise)$shape,edge.color=E(graph_i2.revise)$color2.revise,layout=l)
    
    #plot(graph_i, edge.width=E(graph_i)$width,vertex.frame.width=0.1,vertex.size=V(graph_i)$size, vertex.label=NA,vertex.shapes=V(graph_i)$shape,edge.color=E(graph_i)$color,layout=l)
  
  dev.off()
  
  ###write component vertices and edges data frame into textfile
  
  df.igraph_i<-igraph::as_data_frame(graph_i,what="both")
  write.table(df.igraph_i$vertices, paste0(i,".seed.",seed,".vertices.txt")%>% gsub('\\|','.',.), row.names=FALSE, col.names=TRUE, sep="\t", quote=F)
  #write.table(df.igraph_i$vertices[,c(1:10,16,17,19,20)], paste0(i,".seed.",seed,".vertices.txt")%>% gsub('\\|','.',.), row.names=FALSE, col.names=TRUE, sep="\t", quote=F)
  write.table(df.igraph_i$edges[,1:3], paste0(i,".seed.",seed,".edges.txt")%>% gsub('\\|','.',.), row.names=FALSE, col.names=TRUE, sep="\t", quote=F)
  
  
}
}
######################### The Following Part Is Block Comment #####################

  if (FALSE) {
  file8=paste0(i,".seed.",seed,".LabelPropagation.VS.DirectNeighbor.pdf") %>% gsub('\\|','.',.)
  
  pdf(file=file8,width=18,height=12)
  #pdf(file=file2,width=18,height=12)
  par(mfrow=c(1,2))
  #ceb <- cluster_edge_betweenness(graph_i,weights=NULL)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, graph_i,layout=l,vertex.label=NA,main="Community detection based on edge betweenness")
  
  #cfg <- cluster_fast_greedy(graph_i)
  #plot(cfg, graph_i,layout=l,vertex.label=NA,main="Community detection based on greedy optimization of modularity")
  
  #clp <- cluster_label_prop(graph_i)
  plot(clp, graph_i,layout=l, vertex.label=NA,main="Community detection based on propagating labels")
  
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("Network Order:", length(landscape)),
    paste("Network Size:", ecount(graph_i)),
    paste("Network Membership:", membership),
    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
    paste("Community Enhancers:",community.E ),
    paste("Community Promoters:",community.P ),
    paste("Community Size:",community.size ),
    paste("Community Enhancer Signal:",community.E.signal),
    paste("Community Signal:",community.signal)
  ))
  
  plot(graph_i.p.neighbor2,
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.p.neighbor2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Direct edges")
 
   legend("topleft", legend=c(
    paste("N1 Enhancers:",E ),
    paste("N1 Promoters:",P ),
    paste("N2 Enhancers:",E2 ),
    paste("N2 Promoters:",P2 ),
    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
    paste("N1 Signal:",neighborhood.signal.order1),
    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
    paste("N2 Signal:",neighborhood.signal.order2)
  ))
  
  
  dev.off()
  

  file9=paste0(i,".seed.",seed,".label.propagation.community.preset.neighbors.v2.0",".pdf")
  
  pdf(file=file9,width=18,height=12)
  #pdf(file=file2,width=18,height=12)
  par(mfrow=c(2,2))
  #ceb <- cluster_edge_betweenness(graph_i,weights=NULL)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, graph_i,layout=l,vertex.label=NA,main="Community detection based on edge betweenness")
  
  #cfg <- cluster_fast_greedy(graph_i)
  #plot(cfg, graph_i,layout=l,vertex.label=NA,main="Community detection based on greedy optimization of modularity")
  
  #clp <- cluster_label_prop(graph_i)
  
  plot(graph_i.p.neighbor2,
       layout=l,
       vertex.label.dist=0.6,
       vertex.label=NA,
       edge.color=E(graph_i.p.neighbor2)$edge.color,
       vertex.label.cex=0.7,
       vertex.label.color="black",
       main="Direct edges")
  
  legend("topleft", legend=c(
    paste("Gene:",i),
    paste("Network Order:", length(landscape)),
    paste("Network Size:", ecount(graph_i)),
    paste("Network Membership:", membership),
    paste("Degree:", V(graph_i)$degree[V(graph_i)$name_full==i]),
    paste("N1 Enhancers:",E ),
    paste("N1 Promoters:",P ),
    paste("N2 Enhancers:",E2 ),
    paste("N2 Promoters:",P2 ),
    paste("N1 Network Size:",ecount(graph_i.neighbor1)),
    paste("N2 Network Size:",ecount(graph_i.neighbor2)),
    paste("N1 Enhancer Signal:",neighborhood.E.signal.order1),
    paste("N1 Signal:",neighborhood.signal.order1),
    paste("N2 Enhancer Signal:",neighborhood.E.signal.order2),
    paste("N2 Signal:",neighborhood.signal.order2)
  ))
  
  plot(clp, graph_i,layout=l, vertex.label=NA,main="Default")
  
  legend("topleft", legend=c(
    
    paste("Community Enhancers:",community.E ),
    paste("Community Promoters:",community.P ),
    paste("Community Size:",community.size ),
    paste("Community Enhancer Signal:",community.E.signal),
    paste("Community Signal:",community.signal)
  ))
  
  plot(clp2, graph_i,layout=l, vertex.label=NA,main="Solution 2")
  
  legend("topleft", legend=c(
    
    paste("Community Enhancers:",community.E2 ),
    paste("Community Promoters:",community.P2 ),
    paste("Community Size:",community.size2 ),
    paste("Community Enhancer Signal:",community.E.signal2),
    paste("Community Signal:",community.signal2)
  ))
  
  plot(clp3, graph_i,layout=l, vertex.label=NA,main="Solution 1")
  
  
 
  dev.off()
} 
