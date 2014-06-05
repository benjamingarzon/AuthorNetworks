
rm(list=ls())

library(sna)

workDir <- "/home/benjamingarzon/Data/AuthorNetworks"
source("/home/benjamingarzon/Software/AuthorNetworks/ncbi.search.functions.R")

term <- '\"Quantitative+Susceptibility+Mapping\"'
# retrieve info  
#info <- getInfo(term)
#save(info, file="/home/benjamingarzon/Data/AuthorNetworks/.RNetworkData")
load(paste(workDir, ".RNetworkData", sep="/"))

authorPaperMatrix <- getAuthors(info)
adjacencyMatrix <- authorPaperMatrix %*% t(authorPaperMatrix)
nPapers = diag(adjacencyMatrix)
adjacencyMatrix <- adjacencyMatrix[nPapers>1, nPapers>1]

nPapers = diag(adjacencyMatrix)

colors <- heat.colors(max(nPapers), alpha = 1)
labels <- rownames(adjacencyMatrix)
labels[nPapers < 4] <- seq(length(nPapers))[nPapers < 4]

png(file=paste(workDir, "QSM.png", sep="/"), width = 1500, height = 1500)
gplot(adjacencyMatrix, mode="kamadakawai", arrowhead.cex = 0, displaylabels=TRUE, edge.col="blue", label=labels, label.cex=1.1, gmode="graph", edge.lwd = adjacencyMatrix, vertex.cex=1.5, label.pos=5, vertex.col=colors[nPapers], displayisolates=TRUE, jitter=TRUE)
#gplot(reducedMatrix, mode="fruchtermanreingold", arrowhead.cex = 0, displaylabels=TRUE, edge.col="blue", label.cex=.8, gmode="graph", edge.lwd =.5*reducedMatrix, vertex.cex=2, vertex.col=colors, label.pos=5, displayisolates=FALSE)

dev.off()
authorInfo <- data.frame(Index=seq(length(nPapers)), Author=rownames(adjacencyMatrix), Papers=nPapers)
write.table(authorInfo, file=paste(workDir, "Authors.txt", sep="/"), row.names=FALSE, col.names=TRUE, sep="\t")
