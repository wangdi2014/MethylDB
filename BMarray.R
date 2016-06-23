require(GenomicRanges)
require(bsseq)
#  Create BSseq class for chromosome 1, 2 , ..., 22, X, Y
#  So, I wish we do not use the bsseq structure
#  But, our tcga450K structure constructor Called BMarray

# Class BMarray process Bisulfite Methylation microArray Chip data
# BMarray
# inherited from GenomicsRanges : SummarizedExperiment
setClass("BMarray", contains = "SummarizedExperiment",
         representation(trans="function", parameters="list", annotation = "character"))  # notice, here parameters is used as a public command record

setGeneric("pData", function(object){
  standardGeneric("pData")
})

# Warning，setGeneric declare first and then setMethod
# setGeneric("pData<-", function(object, value){
#   standardGeneric("pData<-")
# })

setGeneric("hasBeensmoothed", function(object){
  standardGeneric("hasBeensmoothed")
})

# Assay selector
setGeneric("selectAssays", function(object, value){
  standardGeneric("selectAssays")
})

# Annotation
setGeneric("getAnnotation", function(object){
  standardGeneric("getAnnotation")
})

setMethod("selectAssays", signature(object="BMarray", value="character"), function(object, value){
  if(!is(try(object@assays, silent = TRUE), "try-error")){
    object@assays$data[[value]];
  }
  else
    stop("error: attribute is not found");
})




# Validation BMarray's method and data obey rules.
setValidity("BMarray", function(object) {
  msg <- validMsg(NULL, .checkAssayNames(object, c("M", "Cov")))
  if(class(rowData(object)) != "GRanges")
    msg <- validMsg(msg, sprintf("object of class %s need to be GRanges in slot rowData", class(object)))
  # Check for Methylation and Coverage ~ Cov = M + U in TCGA lvl-2 data
  if(min(assay(object, "M")) < 0)
    msg <- validMsg(msg, "The M must be positive");
  if(min(assay(object, "Cov")) < 0)
    msg <- validMsg(msg, "The Cov must be positive");
  if(max(assay(object, "M") - assay(object, "Cov")) > 0)
    msg <- validMsg(msg, "The M value must be smaller than Cov value");
  if(is.null(msg)) TRUE else msg
})

# getMethod: colData [biological replicates]
setMethod("pData", "BMarray", function(object){
  object@colData
})

# getGranges rowData
setMethod("granges", signature(x = "GRanges"), function(x){
  x@gr;
})

setMethod("show", "BMarray", function(object){
  cat("An object of type 'BMarray' with \n")
  cat(" ", nrow(object), "methylation loci in 450K array \n")
  cat(" ", ncol(object)/2, "samples with tumor and normal replicates \n")
  if(hasBeensmoothed(object)){
    cat("has been smoothed with \n")
    cat(" ", object@parameters$smoothText, "\n")
  }
  else{
    cat("has not been smoothed . . . \n")
  }
})

# SummarizeExperiment{GenomicsRanges} contains colData of DataFrame
setReplaceMethod("pData", signature(object = "BMarray", 
                                    value = "data.frame"), function(object, value){
                                      colData(object) <- as(value, "DataFrame")
                                    })

# Length of the GRanges object; inherited from SummarizeExperiment
setMethod("length", "BMarray", function(x){
  length(granges(x));
})

# hasBeensmoothed check function
# bool value
setMethod("hasBeensmoothed", "BMarray", function(object){
  "coef" %in% assayNames(object);
})

# Class method : getBeta
# get M, Beta or CN from data
setGeneric("getBeta", function(object, type = c("beta", "M", "CN"), offset = 0, betaThreshold = 0){
  standardGeneric("getBeta")
})

# getBeta For BMarray
setMethod("getBeta", "BMarray", function(object, type = c("beta", "M", "CN"), offset = 0, betaThreshold = 0){
  type = match.arg(type)
  if(type == "beta"){
    Beta =  pmax(getBMarray(object, type = "M") / (getBMarray(object, type = "Cov") + offset), betaThreshold);
    # colnames(Beta) = pData(object)[, "Pair"];
    colnames(Beta) = rownames(pData(object));
    rownames(Beta) = getBMarray(object, type = "gr")@elementMetadata[, "refnames"];
    return(Beta);
  } else if(type == "M"){
    M =  .logit( getBMarray(object, type = "M") / (getBMarray(object, type = "Cov") + offset));
    # colnames(M) = pData(object)[, "Pair"];
    colnames(M) = rownames(pData(object));
    rownames(M) = getBMarray(object, type = "gr")@elementMetadata[, "refnames"];
    return(M);
  }else{
    CN = log2(getBMarray(object, type = "Cov"));
    # colnames(CN) = pData(object)[, "Pair"];
    colnames(CN) = rownames(pData(object));
    rownames(CN) = getBMarray(object, type = "gr")@elementMetadata[, "refnames"];
    return(CN);
  }
})


# getBMarray : get different value from BMarray
# function : inspector function
getBMarray <- function(BMarray, type=c("Cov", "M", "coef", "se.coef", "trans", "parameters", "gr")){
  type = match.arg(type)
  if(type %in% c("M", "Cov"))
    return(assay(BMarray, type))
  if(type %in% c("coef", "se.coef") && type %in% assayNames(BMarray))
    return(assay(BMarray, type))
  if(type %in% c("coef", "se.coef"))
    return(NULL)
  if(type == "trans")
    return(BMarray@trans)
  if(type == "parameters")
    return(BMarray@parameters)
  if(type == "gr")
    return(granges(BMarray))
}

# BMarray constructor
BMarray <- function(Cov = NULL, M=NULL, U=NULL, coef=NULL, se.coef=NULL, trans=NULL, 
                    parameters = NULL, pData = NULL, gr = NULL, pos = NULL, chr = NULL, sampleNames = NULL, element.ref=NULL, annotation=TRUE){
  if(is.null(gr)){
    if(is.null(pos) || is.null(chr) || is.null(element.ref))
      stop("Need pos and chr and cg ref from 450K")
    gr <- GRanges(seqnames=chr, ranges=IRanges(start=pos, width=1), refnames=element.ref)
  }
  
  #	 constructor's correct validation
  if(!is(gr, "GRanges"))
    stop("gr needs to be GRanges")
  if(any(width(gr) != 1))
    stop("CpG site's length need to fix in 1")
  if((is.null(M) && is.null(U))  || (is.null(Cov) && is.null(U)) || (is.null(Cov) && is.null(M)))
    stop("Methylation microarray need (M & Cov) or (M & U) or (U & Cov)")
  
  # Methylation uniform
  if(!is.null(M) && !is.null(Cov)){
    if(!is.matrix(M) || !is.matrix(Cov))
      stop("U & Cov need to be matrices")
    if(length(gr)!=nrow(M) || length(gr)!=nrow(Cov))
      stop("gr and U, Cov should have identity dimensions")
  }
  if(is.null(M)){
    if(!is.matrix(U) || !is.matrix(Cov))
      stop("U & Cov need to be matrices")
    if(length(gr)!=nrow(U) || length(gr)!=nrow(Cov) || ncol(U) != ncol(Cov))
      stop("gr and U, Cov should have identity dimensions")
    M <- Cov - U;
  }
  
  if(is.null(Cov)){
    if(!is.matrix(U) || !is.matrix(M))
      stop("U & M need to be matrices")
    if(length(gr)!=nrow(U) || length(gr)!=nrow(M) || ncol(U) != ncol(M))
      stop("gr and U, M should have identity dimensions")
    Cov <- M + U;
  }
  
  # set gr, M, Cov's names to NULL
  names(gr) <- NULL; names(M) <- NULL; names(Cov) <- NULL;
  
  # pData(colData) should be DataFrame<namespace:IRange>
  if(!is(pData, "DataFrame"))
    pData <- as(pData, "DataFrame")
  
  # change data into continous by chromosome --> reduce from GenomicsRanges
  grR <- reduce(gr, with.revmap = TRUE)
  
  # reduce will reduce gr's length, but we wish to recover it to its original length
  mm <- unlist(elementMetadata(grR)[[1]])  # cg sequence list
  gr <- gr[mm]
  M <- M[mm,,drop = FALSE]; Cov <- Cov[mm,,drop = FALSE]
  if(!is.null(coef))
    coef <- coef[mm,,drop = FALSE]
  if(!is.null(se.coef))
    se.coef <- se.coef[mm,,drop = FALSE]
  
  # sampleNames use to M data and Cov data
  if(length(sampleNames) != ncol(M)){
    stop("sample's lenght need to identical to M matrix and Cov matrix's dimensions");
  }
  colnames(M) <- sampleNames
  colnames(Cov) <- sampleNames
  
  if(!is.null(coef)){
    if(ncol(coef) != ncol(M))
      stop("coef need to identical to M matrix's dimension")
    colnames(coef) <- sampleNames
  }
  
  if(!is.null(se.coef)){
    if(ncol(se.coef) != ncol(M))
      stop("se.coef need to identical to M matrix's dimension")
    colnames(se.coef) <- sampleNames
  }
  
  if(length(sampleNames) == nrow(pData))
    rownames(pData) <- sampleNames
  
  assays <- SimpleList(M = M, Cov = Cov, coef = coef, se.coef = se.coef)
  # delete null item in SimpleList; Notice SimpleList only in SummarizExperiment's data of Assays
  assays <- assays[!sapply(assays, is.null)]
  
  # pData format should be represented as the C_n, N_n Type Pair
  if(is.null(pData) || all(dim(pData) == c(0, 0)))
    BMarray <- SummarizedExperiment(rowData = gr, assays = assays)
  else
    BMarray <- SummarizedExperiment(rowData = gr, colData = pData, assays = assays)
  
  BMarray <- as(BMarray, "BMarray")  # binding methods declared above : down-binding class from superclass
  if(is.function(trans))
    BMarray@trans = trans
  if(is.list(parameters))
    BMarray@parameters = parameters
  if(annotation)
    BMarray@annotation = unlist(list(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))
  BMarray
}

# update from BMarray to BSseq
# This is just for test ...
setMethod("updateObject", "BMarray", function(object, ...){
  if(!is(try(object@assays, silent = TRUE), "try-error")){
    return(BSseq(M=selectAssays(object, "M"), Cov=selectAssays(object, "Cov"), pData = pData(object), gr = granges(object)))
  }
  else{
    stop("object did not match")
  }
})

# get cgXXXXX's annotation information
# Args : BMarray object
setMethod("getAnnotation", "BMarray", function(object){
  object <- get(.getAnnotation(object));
  object
})


.checkAssayNames <- function(object, names){ 
  nms <- assayNames(object) 
  if(!all(names %in% nms)) 
    return(sprintf("object of class '%s' needs to have assay slots with names '%s'",class(object), paste0(names, collapse = ", "))) 
  else 
    NULL 
} 

# .getAnnotation
.getAnnotation <- function(object){
  if(all(c("array", "annotation") %in% names(object@annotation)))
    return(sprintf("%sanno.%s", object@annotation["array"], object@annotation["annotation"]))
  stop("unable to get annotation database")
}

# utils functions: will be put into util.R later.
# assayName : return colnames of data collected in data slot
# .logit : to calculate log2(M/(1-M))
assayNames <- function(object){ names(object@assays$field("data")) }

.logit <- function(x, offset = 100){ log2(x / (1 - x)) }

.unittest <- function(){
  # BMarray <- function(Cov = NULL, M=NULL, U=NULL, coef=NULL, se.coef=NULL, trans=NULL, 
  #        parameters = NULL, pData = NULL, gr = NULL, pos = NULL, chr = NULL, sampleNames = NULL, element.ref=NULL)
  M <- matrix(0:8, 3, 3)
  Cov <- matrix(1:9, 3, 3)
  pdata <- data.frame(Type = c("cancer", "cancer", "cancer"), Pair = c("pair1", "pair2", "pair3"))
  sampleNames <- c("A1", "A2", "A3")
  chr <- c("chr1", "chr2", "chr3")
  pos <- 1:3
  element <- c("cg1", "cg2", "cg3")
  BMtmp <- BMarray(Cov = Cov, M = M, pData = pdata, chr = chr, pos = pos, sampleNames = sampleNames, element.ref = element)
}
