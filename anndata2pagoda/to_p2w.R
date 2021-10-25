suppressMessages(library(pagoda2))
suppressMessages(library(Matrix))

to_p2w <- function(adata,use_rep,pals,title,output){
    
    additionalMetadata <- list()
    ## for Infomap use hue values from 0.1 to 0.5
    
    for (p in names(pals)){
        additionalMetadata[[p]] <- p2.metadata.from.factor(adata@colData[[p]], 
                                                                   displayname = p,
                                                                   pal=pals[[p]])
    }
    
    cnts=adata@assays@data$X
    rownames(cnts)=rownames(adata)
    colnames(cnts)=colnames(adata)
    cat("    generating pagoda2 object\n")
    p2=Pagoda2$new(cnts,verbose=FALSE)
    p2$counts=t(cnts)
    p2$depth=adata@colData$total_counts

    p2$adjustVariance(verbose=FALSE)
    p2$reductions$PCA=reducedDim(adata)
    p2$clusters$PCA$leiden=adata@colData$leiden
    p2$embeddings$PCA$UMAP=reducedDim(adata,use_rep)
    p2$n.cores=1
    cat("    getting Hierarchical Diff Expression Aspects\n")
    
    leiden="leiden"
    if ("leiden_noname" %in% names(adata@colData)){
        leiden="leiden_noname"
        p2$clusters$PCA[[leiden]]=adata@colData$leiden_noname
    }
    hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA', clusterName=leiden,verbose=F)
    
    suppressMessages(library(org.Mm.eg.db))
    ids <- unlist(lapply(mget(colnames(p2$counts), org.Mm.egALIAS2EG, ifnotfound=NA), function(x) x[1]))
    # reverse map
    rids <- names(ids)
    names(rids) <- ids
    # list all the ids per GO category
    go.env <- list2env(eapply(org.Mm.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
    genesets <- hierDiffToGenesets(hdea)
    library(GO.db)
    termDescriptions <- Term(GOTERM[names(go.env)]) # saves a good minute or so compared to individual lookups
    sn <- function(x) { names(x) <- x; x}  ## utility function

    genesets.go <- lapply(sn(names(go.env)),function(x) {
      list(properties=list(locked=TRUE, genesetname=x, shortdescription=as.character(termDescriptions[x])), genes=c(go.env[[x]]))
    })

    ## concatenate
    genesets <- c(genesets, genesets.go)
    deSets <- get.de.geneset(p2, groups = p2$clusters$PCA[[leiden]], prefix = 'de_')
    ## concatenate
    genesets <- c(genesets, deSets)
    appmetadata <- list(apptitle = title)
    cat("    making Gene Knn graph\n")
    p2$makeGeneKnnGraph(n.cores = 1,verbose=FALSE)

    ## # Make a list for our metadata
    

    p2w <- make.p2.app(p2, 
        dendrogramCellGroups = p2$clusters$PCA[[leiden]],
        additionalMetadata = additionalMetadata,
        geneSets = genesets,
        appmetadata = appmetadata,
        show.clusters = FALSE # Hide the clusters that were used for the dendrogram from the metadata
      )
    cat("    saving p2 web bin file\n")
    p2w$serializeToStaticFast(output)
}