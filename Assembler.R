# This script is for parallel processing and create pairwised or non-pairwised 
# BMarray data.

load("../data/methyl450patient.merge.rda")  # Data for methyl450kpatient
                                            # CpGref  --- CpG reference across all patients
                                            # Patient --- All patients' ID, suffix of rda
                                            # Pair    --- All paired patients' ID
library(parallel)
source("./BMarray.R")
## ifndef default savedirectory
savedirect.default = "../data";

## function    : Assembler
## description : Build BMarray object for specific cancer's database
## @arguments  : db450kDirectory --- 450k.process.R's output : .rda
##               mcore           --- multi-core process here
##               paired          --- TRUE for pairwised searching
##                                   FALSE for totally searching
## @output     : BMarray data structure
AssembleBMarray <- function(db450kDirectory = NULL, mcore = 4, paired = FALSE) {
  options(warn = -1)
  if(is.null(db450kDirectory) || mcore < 1)
    stop("need db450kDirectory and how many cores do you wish to use?")
  # load patient datafram
  load(paste(db450kDirectory, "methyl450patient.merge.rda", sep="/"))   
  # ==> methyl450kpatient list(CpGref=CpG.ref, Patient=Patient.ref, Pair=Pair.ref)
  # load 450k metadata
  load(paste(db450kDirectory, "BRCA_metadb.rda", sep="/"))
  # tabulate M/Cov/CN/Sample Name
  # initialization
  CpG.ref <- methyl450kpatient$CpGref
  Patient.ref <- NULL; Patient.size <- NULL;
  if(paired == TRUE){
    cat(sprintf("89 paired samples are assembled ...\n"))
    Patient.ref <- methyl450kpatient$Pair
    Patient.size <- length(Patient.ref)
  }
  else{
    cat(sprintf("all samples are assembled ... \n"))
    Patient.ref <- methyl450kpatient$Patient
    Patient.size <- length(Patient.ref)
  }
  # multi-core preparation
  mc <- getOption("mc.cores", mcore)
  ############################  filter process 
  
  cat(sprintf("calculate tumor ... \n"))
  ## tumor result -- For non-pairwise; there maybe NULL in result.tumor
  if(paired == TRUE){
    result.tumor <- mclapply(Patient.ref, function(x){
      m <- Methyl_batchretrieve(type = "tumor", sampleidx = x, cpglist = CpG.ref, db450kDirectory = db450kDirectory, keep_CN = TRUE); 
    }, mc.cores = mc);
  }
  else{
    result.tumor <- lapply(Patient.ref, function(x){
      m <- Methyl_batchretrieve(type = "tumor", sampleidx = x, cpglist = CpG.ref, db450kDirectory = db450kDirectory, keep_CN = TRUE); 
    });
  }
  names(result.tumor) <- Patient.ref
  ## remove NULL
  result.tumor[sapply(result.tumor, is.null)] <- NULL
  patient.idx = names(result.tumor)
  ## Build Cov M matrix for tumor list
  if(paired == TRUE){
    CpGidx = Reduce(intersect, mclapply(result.tumor, function(x){ x$Element.Ref }, mc.cores = mc));
  }
  else{
    CpGidx = Reduce(intersect, lapply(result.tumor, function(x){ x$Element.Ref }));
    mcore = 1;
  }
  M.tumor = .mapreduce(oplst = result.tumor, mcore = mcore, map = "Element.Ref", rules = CpGidx, by = "M", type = "col")
  Cov.tumor = .mapreduce(oplst = result.tumor, mcore = mcore, map = "Element.Ref", rules = CpGidx, by = "Cov", type = "col")
  SampleName.tumor =  paste(patient.idx, "T", sep="."); 
  element = CpGidx;
  # data tmp kept
  tumor.length = length(result.tumor)
  tumor.name = names(result.tumor)
  rm(result.tumor)
  cat(sprintf("calculate normal ... \n"))
  ## normal resulr -- For non-pairwise; there maybe NULL in result.normal
  if(paired == TRUE){
    result.normal <- mclapply(Patient.ref, function(x){
      m <- Methyl_batchretrieve(type = "normal", sampleidx = x, cpglist = CpG.ref, db450kDirectory = db450kDirectory, keep_CN = TRUE); 
    }, mc.cores = mc);
  }
  else{
    result.normal <- lapply(Patient.ref, function(x){
      m <- Methyl_batchretrieve(type = "normal", sampleidx = x, cpglist = CpG.ref, db450kDirectory = db450kDirectory, keep_CN = TRUE); 
    });
  }
  names(result.normal) <- Patient.ref
  ## remove NULL
  result.normal[sapply(result.normal, is.null)] <- NULL
  patient.idx = names(result.normal)
  ## Build Cov M matrix for normal list
  if(paired == TRUE){
    CpGidx = Reduce(intersect, mclapply(result.normal, function(x){ x$Element.Ref }, mc.cores = mc));
  }
  else{
    CpGidx = Reduce(intersect, lapply(result.normal, function(x){ x$Element.Ref }));
    mcore = 1;
  }
  M.normal = .mapreduce(oplst = result.normal, mcore = mcore, map = "Element.Ref", rules = CpGidx, by = "M", type = "col")
  Cov.normal = .mapreduce(oplst = result.normal, mcore = mcore, map = "Element.Ref", rules = CpGidx, by = "Cov", type = "col")
  SampleName.normal = paste(patient.idx, "N", sep=".")
  
  ## Create merged M and Cov and pData for BMarray object
  M = do.call(cbind, list(M.tumor, M.normal))
  Cov = do.call(cbind, list(Cov.tumor, Cov.normal))
  sampleNames = c(SampleName.tumor, SampleName.normal)
  chr = as.character(result.normal[[1]]$Chr)
  pos = result.normal[[1]]$start;
  
  ## pData created
  #pData <- data.frame(Type = c(rep("tumor", length(result.tumor)), rep("normal", length(result.normal))), 
  #                    Pair = c(names(result.tumor), names(result.normal)))
  pData <- data.frame(Type = c(rep("tumor", tumor.length), rep("normal", length(result.normal))), 
                      Pair = c(tumor.name, names(result.normal)))
  bmarray.result = BMarray(Cov = Cov, M = M, pData = pData, chr = chr, pos = pos, sampleNames = sampleNames, element.ref = CpGidx)
  
  return(bmarray.result)
}

## function    : Methyl_batchretrieve
## description : batch retrieval of illunima 450k Rda data for matrix generation for BMarray
## @arguments  : type            -- which cancer type select : tumor/normal
##               sampleidx       -- patiente based 450k methylation database's id : sampleidx.rda
##               cpglist         -- the CpG island list which is uniform across all sample's database
##               db450kDirectory -- 450k methylation database's directory
##               keep_CN         -- copy number keep or not [bool : TRUE || FALSE]
## @output     : list (M, Cov, start)
Methyl_batchretrieve <- function(type = c("normal", "tumor"), 
                                 sampleidx = NULL, cpglist = NULL, 
                                 db450kDirectory = NULL, keep_CN = FALSE) {
  options(warn = -1)
  if(is.null(sampleidx) || is.null(cpglist) || is.null(db450kDirectory))
    stop("ampleID and directory information should be proffered sametime")
  type = match.arg(type)
  cat(sprintf("%s \n", sampleidx));
  # load patient specific database : 
  # .rda ==> merge.methyl ~ list(tumor = tumor_lst, paired = paired, names = tumor_lst[,1])
  load(paste(db450kDirectory, paste(sampleidx, "rda", sep="."), sep = "/"))
  if(type %in% names(merge.methyl)){
    methyl.df <- merge.methyl[[type]];
    rownames(methyl.df) <- as.character(methyl.df$Composite.Element.REF);
    # reorder by cpglist
    methyl.df = na.omit(methyl.df[cpglist,])
    rownames(methyl.df) <- NULL;
    result = NULL
    rm(merge.methyl)
    # keep_CN or not
    if(keep_CN){
      result = data.frame(Element.Ref = as.character(methyl.df$Composite.Element.REF),
                          M = methyl.df$Methylated_Intensity,
                          Cov = methyl.df$Methylated_Intensity+methyl.df$Unmethylated_Intensity,
                          CN = log2(methyl.df$Methylated_Intensity+methyl.df$Unmethylated_Intensity),
                          Chr = paste("chr", as.character(methyl.df$Chromosome),sep=""),
                          start = methyl.df$Genomic_Coordinate);
    }
    else{
      result = data.frame(Element.Ref = as.character(methyl.df$Composite.Element.REF),
                          M = methyl.df$Methylated_Intensity,
                          Cov = methyl.df$Methylated_Intensity+methyl.df$Unmethylated_Intensity,
                          Chr = paste("chr", as.character(methyl.df$Chromosome),sep=""),
                          start = methyl.df$Genomic_Coordinate);
    }
  }
  else{
    result = NULL;
  }
  return(result);
}

## function    : .mapreduce
## description : a functioanl perspective to filter result of a List contain 450k data
## @arguments  : oplst  --- operated list
##               mcore  --- mcore in parallel package
##               map    --- filter targets (what will be reduced by rules)
##               rules  --- filter rules [based on maps] 
##               by     --- column will be kept
##               type   --- result orginized by Col | Row
## @output     : a matrix
.mapreduce <- function(oplst = NULL, mcore = mcore, map = NULL, rules = NULL, by = NULL, type = c("col", "row")) {
  options(warn = -1)
  if(is.null(oplst) || is.null(map) || is.null(rules) || is.null(by))
    stop("map & rules & oplst & by are needed")
  type = match.arg(type)
  if(mcore >= 2){
    mc = getOption("mc.cores", mcore)
    if(type == "col"){
      result <- do.call(cbind, mclapply(oplst, function(x){
                                                rownames(x) = x[[map]]
                                                x <- x[rules, ];    rownames(x) <- NULL
                                                return(x[[by]])
                                               }, mc.cores = mc));
    }
    else if(type == "row"){
      result <- do.call(rbind, mclapply(oplst, function(x){
                                                rownames(x) = x[[map]]
                                                x <- x[rules, ];    rownames(x) <- NULL
                                                return(x[[by]])
                                               }, mc.cores = mc));
    }
    else{
      stop("operation is not allowed")
    }
  }
  else{
    if(type == "col"){
      result <- do.call(cbind, lapply(oplst, function(x){
        rownames(x) = x[[map]]
        x <- x[rules, ];    rownames(x) <- NULL
        return(x[[by]])
      }));
    }
    else if(type == "row"){
      result <- do.call(rbind, lapply(oplst, function(x){
        rownames(x) = x[[map]]
        x <- x[rules, ];    rownames(x) <- NULL
        return(x[[by]])
      }));
    }
    else{
      stop("operation is not allowed")
    }
  }
  return(result)
}

.unittest_Assembler <- function(){
  result.paired <- AssembleBMarray(db450kDirectory = "../data", mcore = 6, paired = TRUE)
  result.unpaired <- AssembleBMarray(db450kDirectory = "../data", mcore = 6, paired = FALSE)
  CN.paired = log2(result.paired@assays$data$Cov);
  CN.unpaired = log2(result.unpaired@assays$data$Cov);
  tmp <- list(BMarray = result.paired, CN = CN.paired)
  save(tmp, file = "../data/BMarrayDB.rda")
  tmp <- list(BMarray = result.unpaired, CN = CN.unpaired)
  save(tmp, file = "../data/BMarrayDB_all.rda")
}







