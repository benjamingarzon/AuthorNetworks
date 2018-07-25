
rm(list=ls())

library(sna)
library(easyPubMed)
workDir <- "~/projects/AuthorNetworks"

# retrieve info  
my_query <- '\"delay discounting\"[All Fields] OR \"intertemporal choice\"[All Fields] OR \"temporal discounting\"[All Fields] OR \"time preference\"[All Fields]'
batches <- batch_pubmed_download(pubmed_query_string = my_query, dest_file_prefix = "dd_sample")
# Retrieve the full name of the XML file downloaded in the previous step

df = NULL
for (i in 1:length(batches))
 df <- rbind(df, table_articles_byAuth(pubmed_data = batches[i], max_chars = 0))
             
# Alternatively, the output of a fetch_pubmed_data() could have been used
#
# Printing a sample of the resulting data frame
df$address <- substr(df$address, 1, 28)
df$jabbrv <- substr(df$jabbrv, 1, 9)
df$name <- paste(df$firstname, df$lastname)
print(df[1:10, c("pmid", "year", "jabbrv", "name")])  

authors <- unique(df$name)    
PMID = df$pmid
unPMID <- unique(PMID)
authorPaperMatrix <- matrix(data=0, nrow=length(authors), ncol=length(unPMID), dimnames=list(authors, unPMID))

for (i in (1:length(authors)))
  authorPaperMatrix[authors[i], PMID[df$name==authors[i]]] <- 1 

authorPaperMatrix <- authorPaperMatrix[rownames(authorPaperMatrix)!="NA NA", ]

minPapers <- 10
adjacencyMatrix <- authorPaperMatrix %*% t(authorPaperMatrix)
nPapers = diag(adjacencyMatrix)
adjacencyMatrix <- adjacencyMatrix[nPapers>minPapers, nPapers>minPapers]

nPapers = diag(adjacencyMatrix)

colors <- heat.colors(max(nPapers), alpha = 1)
labels <- rownames(adjacencyMatrix)
#labels[nPapers < 4] <- seq(length(nPapers))[nPapers < 4]

png(file=paste(workDir, "DelayDiscounting.png", sep="/"), width = 1500, height = 1500)
gplot(adjacencyMatrix, mode="fruchtermanreingold", arrowhead.cex = 0, 
      displaylabels=TRUE, edge.col="blue", label=labels, label.cex=1.1, 
      gmode="graph", edge.lwd = adjacencyMatrix, vertex.cex=1.5, 
      label.pos=5, vertex.col=colors[nPapers], displayisolates=TRUE, jitter=TRUE)
#gplot(reducedMatrix, mode="fruchtermanreingold", arrowhead.cex = 0, displaylabels=TRUE, edge.col="blue", label.cex=.8, gmode="graph", edge.lwd =.5*reducedMatrix, vertex.cex=2, vertex.col=colors, label.pos=5, displayisolates=FALSE)

dev.off()
authorInfo <- data.frame(Index=seq(length(nPapers)), Author=rownames(adjacencyMatrix), Papers=nPapers)
write.table(authorInfo, file=paste(workDir, "Authors.txt", sep="/"), row.names=FALSE, col.names=TRUE, sep="\t")
