library(compiler)
library(scrapeR)
library(RNCBI)
library(sna)

# search for a term in a database, return ID by publication date
doSearch <- function(term, db="pubmed", retstart=0, retmax=500){
  ncbi<-NCBI()
  
  esearch <- ESearch(ncbi)
  esearch <- setRequestParameter(esearch, "db", db)
  esearch <- setRequestParameter(esearch, "term", term)
  esearch <- setRequestParameter(esearch, "sort", "PublicationDate")
  esearch <- setRequestParameter(esearch, "retmax", retmax)
  esearch <- setRequestParameter(esearch, "retstart", retstart)
  esearch <- requestESearch(esearch)
  esearchRes <- getResults(esearch)
  
  IDlist <- unlist(esearchRes$idlist$id)
}

noMidName <- function(name){
  names <- strsplit(name, " ")[[1]]
  return(paste(names[1], names[length(names)], sep=" "))
  
}

fixScandinavian <- function(page){  
  page <- gsub("Ã¶","ö",page)
  page <- gsub("Ã¤","ä",page)
  page <- gsub("Ã¢","å",page)  
}

# return the first descriptor
returnFirst <- function(l){
  return(xmlValue(l[[1]]))
  
}

# returns a list with authors, Publication date, ...
getInfo <- function(term, db="pubmed", maxResultsPerSearch=500){
  
  info <- list()
  j <- 1
  jj <- 1
  offset <- 0 
  retmax <- nResults <- maxResultsPerSearch
  
  while(nResults == retmax){ 
    
    
    IDlist <- doSearch(term, retstart=offset, retmax=retmax)
    offset <- offset + retmax
    nResults <- length(IDlist)
    url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=", db, "&retstart=", 0, "&retmax=", retmax,"&retmode=xml&id=", paste(IDlist,collapse=","), sep="")
    page <- getURL(url, .encoding="UTF-8")    
    doc <- xmlTreeParse(page, useInternal = TRUE, encoding="UTF-8")
    
    top <- xmlRoot(doc)
    
    # can be done more flexible!!
    
    for (i in 1:length(names(top))){
      pmid <- xmlValue(top[[i]][["MedlineCitation"]][["PMID"]])    
      
      history <- top[[i]][["PubmedData"]][["History"]]
      if(is.null(history)) next 
      
      articleDate <- top[[i]][["PubmedData"]][["History"]][[which(xmlSApply(history,xmlGetAttr,"PubStatus")=="pubmed")]]
      authorList <- top[[i]][["MedlineCitation"]][["Article"]][["AuthorList"]]
      articleTitle <- top[[i]][["MedlineCitation"]][["Article"]][["ArticleTitle"]]
      meshHeadingList <- top[[i]][["MedlineCitation"]][["MeshHeadingList"]]
      
      if(!is.null(meshHeadingList)){
        
        keywords <- unlist(xmlApply(meshHeadingList, returnFirst))
        names(keywords) <- NULL
      } else {
        keywords <- NULL
      }
      
      
      # removing books, for example
      if(is.na(pmid) | is.null(authorList) | is.null(articleDate)){
        jj <- jj+1
        next
      }
      
      info[[j]] <- list(ArticleDate=xmlSApply(articleDate, xmlValue),
                        AuthorList=xmlSApply(authorList, xmlSApply, xmlValue),
                        ArticleTitle=xmlValue(articleTitle), PMID=pmid, keywords=keywords )
      
      j <- j+1  
    }
  }
  
  cat("Total: ", length(info), "\n")
  
  return (info)  
}

# given a list of papers, return a list of corresponding TACs
calculateTAC <- function(info ){
  l <- length(info)
  year <- rep(0,lenght=l)
  for (i in (1:l)){
    print(i)
    year[i] <- as.numeric(info[[i]]$ArticleDate["Year"])
    
  }
  counts <- table(year)
  year <- as.numeric(names(counts))
  return(list(year=year, counts=counts))
}

# plot a list of TACs
plotTAC <- function(TAC, color="red"){
  lines(TAC$year,TAC$counts,col=color)    
}

# get all the authors
getAuthors <- function(info){
  
  l <- length(info)
  name <- PMID <- index <- unPMID <- NULL
#  unPMID <- rep("", length=l)
  k <- 1
  ii <- 1
  for (i in (1:l)){
    
    authors <- info[[i]]$AuthorList
    
    if (typeof(authors)=="character") {
      if (is.null(dim(authors))) next
      
      unPMID <- c(unPMID, info[[i]]$PMID)
      ii <- ii + 1
      for (j in 1:dim(authors)[2]){
        #  cat(i,j,k)
        name <- c(name, noMidName(fixScandinavian(paste(authors[1,j], authors[2,j], sep=" "))))
        PMID <- c(PMID, info[[i]]$PMID)
        
        index <- c(index, i )
        k <- k+1
      }
    }
    else{
      if (is.null(length(authors))) next
      
      unPMID <- c(unPMID, info[[i]]$PMID)
      ii <- ii + 1
      for (j in 1:length(authors)){
        #  cat(i,j,k)
        name <- c(name, noMidName(fixScandinavian(paste(authors[j]$Author["ForeName"], authors[j]$Author["LastName"], sep=" "))))
        PMID <- c(PMID, info[[i]]$PMID)
        
        index <- c(index, i )
        k <- k+1
      }
      
    }

  }
  
  # possibly needs identification of authors with different spelling
  authors <- unique(name)    
  
  authorPaperMatrix <- matrix(data=0, nrow=length(authors), ncol=length(unPMID), dimnames=list(authors, unPMID))
  # parallelize !!
  for (i in (1:length(authors))){
    #print(i)
    PMID[name==authors[i]]
    authorPaperMatrix[authors[i], PMID[name==authors[i]]] <- 1 
  }
  authorPaperMatrix <- authorPaperMatrix[rownames(authorPaperMatrix)!="NA NA", ]
  return(authorPaperMatrix)
}


getAuthorsCmp <- cmpfun(getAuthors)

calculateParameters <- function(authorPaperMatrix, adjacencyMatrix){
  
  nPapers <- ncol(authorPaperMatrix)
  nAuthors <- nrow(authorPaperMatrix)
  
  meanDegree <- mean(degree(adjacencyMatrix))
  
  nComponents <- components(adjacencyMatrix)
  central <- centralization(adjacencyMatrix, degree)
  
  authorsPerComponent <- nAuthors/nComponents
  
  result <- c(nPapers = nPapers, nAuthors = nAuthors, meanDegree = meanDegree, nComponents = nComponents, centralization = central, authorsPerComponent = authorsPerComponent)
  return(result)
}


calculateParametersCmp <- cmpfun(calculateParameters)

getKeywords <- function(term){
  
  info <- getInfo(term)
  keywords <- NULL
  for (i in 1:length(info)){    
    keywords <- c(keywords, info[[i]]$keywords)
  }
  
  return(paste(keywords, collapse = " / "))  
}