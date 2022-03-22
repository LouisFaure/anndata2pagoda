suppressMessages(library(pagoda2))
suppressMessages(library(Matrix))

to_p2w <- function(adata,use_rep,clustering,pals,title,output){
    
    
    
    cnts=adata@assays@data$X
    rownames(cnts)=rownames(adata)
    colnames(cnts)=colnames(adata)
    cat("    generating pagoda2 object\n")
    p2=Pagoda2$new(cnts,verbose=FALSE)
    p2$counts=t(cnts)
    p2$depth=adata@colData$total_counts

    p2$adjustVariance(verbose=FALSE)
    p2$reductions$PCA=reducedDim(adata)
    p2$embeddings$PCA$UMAP=reducedDim(adata,use_rep)
    
    additionalMetadata <- list()
    ## for Infomap use hue values from 0.1 to 0.5
    
    for (p in names(pals)){
        additionalMetadata[[p]] <- p2.metadata.from.factor(adata@colData[[p]], 
                                                                   displayname = p,
                                                                   pal=pals[[p]])
        p2$clusters$PCA[[p]]=adata@colData[[p]]
    }
    
    p2$n.cores=1
    
    clusters <- clustering
    if (paste0(clusters,"_noname") %in% names(adata@colData)){
        clusters <- paste0(clusters,"_noname")
    }
    p2$clusters$PCA[[clustering]] <- adata@colData[[clusters]]
    
    cat("    getting Hierarchical Diff Expression Aspects\n")
    hdea <- p2$getHierarchicalDiffExpressionAspects(type='PCA', clusterName=clustering,verbose=F)
   
    # list all the ids per GO category
    genesets <- hierDiffToGenesets(hdea)
    deSets <- get.de.geneset(p2, groups = p2$clusters$PCA[[clustering]], prefix = 'de_')
    ## concatenate
    genesets <- c(genesets, deSets)
    appmetadata <- list(apptitle = title)
    cat("    making Gene Knn graph\n")
    p2$makeGeneKnnGraph(n.cores = 1,verbose=FALSE)

    ## # Make a list for our metadata
    

    p2w <- make.p2.app(p2, 
        dendrogramCellGroups = p2$clusters$PCA[[clustering]],
        additionalMetadata = additionalMetadata,
        geneSets = genesets,
        appmetadata = appmetadata,
        show.clusters = FALSE # Hide the clusters that were used for the dendrogram from the metadata
      )
                              
    ext=stringr::str_split(output,"[.]")[[1]][2]
    if (ext=="bin"){
        cat("    saving p2 web bin file\n")
        p2w$serializeToStaticFast(output)
    } else if (ext=="RData"){
        cat("    saving p2 web RData file\n")
        save(p2w,file=output)
    }
    
}