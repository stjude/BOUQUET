
 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
 calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
# 	if(drawPlot){  #if TRUE, draw the plot
# 		plot(1:length(inputVector), inputVector,type="l",...)
# 		b <- y_cutoff-(slope* xPt)
# 		abline(v= xPt,h= y_cutoff,lty=2,col=8)
# 		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
# 		abline(coef=c(b,slope),col=2)
# 		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
# 		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
# 	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}

########## PLOT 3D SE SWOOSH


args <- commandArgs()

print('THESE ARE THE ARGUMENTS')
print(args)

#ARGS
outPrefix = args[3]
SignalFile = args[4]
highlightGenes= args[5]

# R --no-save output_PriMROSE RegionsCollapsed.cis.txt RegionsCollapsed.bed.MED1 RegionsCollapsed.bed.Input < PriMROSE_callSuper.R


SignalTable<-read.table(SignalFile, header=F, sep="\t")

SignalTable[SignalTable[,2]<0,2]<-0
cutoff_options<-calculate_cutoff(SignalTable[,2])
SignalTable<-SignalTable[order(-SignalTable[,2]),]
SignalTable$Rank<-rank(SignalTable[,2])
SignalTable$isSuper<-SignalTable[,2] > cutoff_options$absolute


pdf(paste(outPrefix, "_swoosh_", Sys.Date(), ".pdf", sep=""))
plot(SignalTable[order(SignalTable[,2]),2], type='l', col="red", ylab="Dedicated Enhancer Signal", xlab="Gene Ranking")
legend("topleft", legend=c(
  paste("Cutoff:", cutoff_options$absolute),
  paste("Supers:", sum(SignalTable$isSuper)),
  paste("Typicals:", sum(!SignalTable$isSuper))
))
abline(h=cutoff_options$absolute, lty=2)


if(!is.na(highlightGenes)) {
  highlightGenes<-t(do.call('rbind', strsplit(as.character(highlightGenes),',',fixed=TRUE)))
  cols<-rainbow(length(highlightGenes))

  for( i in 1:length(highlightGenes)) {  
    geneID=unlist(strsplit(highlightGenes[i],'|',fixed=TRUE))[1]
    txID=unlist(strsplit(highlightGenes[i],'|',fixed=TRUE))[2]
    alternateID<-paste0(txID,"|",geneID)
    subset<-SignalTable[(SignalTable[,1] == highlightGenes[i])|(SignalTable[,1] == alternateID),]
    print(highlightGenes[i])
    print(subset)
    points(subset$Rank, subset[,2], col=cols[i], pch=19, cex=1.3)

    text(subset$Rank-3000, subset[,2], labels=subset[,1], cex=1, pos=3)
  }
  legend("left", legend=highlightGenes, col=cols, pch=19)

 } 

colnames(SignalTable)<-c("GENE_ID", "ENHANCER_SIGNAL", "GENE_RANK", "IS_SUPER")
write.table(SignalTable, paste(outPrefix, ".AllGenes.table.txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t", quote=F)
write.table(SignalTable[SignalTable$IS_SUPER,], paste(outPrefix, ".SuperGenes.table.txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t", quote=F)

dev.off()
