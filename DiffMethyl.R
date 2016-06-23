# this script is for two variables
# DMP -- Differential Methylation Position
# DMR -- Differential Methylation Region
# DMP's search : glm & F-test  ~ DMP_finder
# DMR's search : bumphunter    ~ DMR_finder
require(limma)
require(stringr)
require(siggenes)
## For unit test, and our analysis, I must introduce some ultility functions for BMarray structure's analysis
source("./BMarray.R")
source("./Utils.R")        ## Utility functions
source("./Annotation.R")   ## Annotation functions 
load("data/BRCA_metadb.rda")


## function    : DMP_engine
## description : find out differential methylation positions by linear regression
##               Y_{ij} = \mu(t_j) + \beta(t_j)X_i + \sum_{k=1}^K \gamma(t_j) Z_{i,k} + \epsilon_{ij}
##               i is for individule ~ sample
##               t_j is for CpG site
##               @reference : Irizarry, Rafael A., et al. "The human colon cancer methylome shows similar hypo-and hypermethylation at conserved tissue-specific CpG island shores." Nature genetics 41.2 (2009): 178-186.
##               based on this paper, we set p-value cutoff to 0.05, and fdr to 0.01, 0.05, 0.10
## @arguments  : dat      --- matrix of M-value or Beta-value  [cpg_sites x sample_size]
##               pheno    --- variable to be observed [grouped factor or continous variable]
##               cofactor --- other factors which will affect output of linear regression
##                            e.g. in Methylation data; age is highly correlated with methylation status
##               type     --- pheno data's type : "categorical" or "continuous"
##               FDR      --- false discovery rate; we should control FDR within a threshold
##               pval     --- significance of positions's methylation status
## @output     : Dmp.tab  --- Differential methylation positions table
DMP_engine <- function(dat, pheno = NULL, cofactor = NULL, type = c("categorical", "continuous"), FDR = 0.05, pval = 0.05) {
  options(warn = -1)
  type = match.arg(type)
  if(is.null(pheno))
    stop("DMP_engine needs pheno ... ")
  if(!is.matrix(dat))
    stop("DMP_engine needs matrix-like data to transfer each-cpg's methylation measurement")
  if(!type %in% c("categorical", "continuous"))
    stop("pheno's type must be constrained in (categorical, continuous)")
  if(is.factor(cofactor))
    cofactor <- .factor2num(cofactor)
  if(!is.vector(cofactor))   # if cofactor is not a vector but a matrix
    stop("DMP engine does not support multiple co-factors now")
  # n return sample size: n_1 tumor samples and n_2 normal samples; n_1 + n_2 = n
  n = ncol(dat)
  r = nrow(dat)
  if(type == "categorical" && !is.null(cofactor)) {    # age or gender cofactor should be considered 
    cat(sprintf("cofactor is set as continuous variables ...\n"))
    pheno = as.factor(pheno);
    dd = data.frame(pheno = pheno, cofactor = cofactor)
    cofactor.design = model.matrix(~ cofactor, dd)
    fit <- lmFit(dat, cofactor.design)
    # leave cofactor-related part out
    beta = fit$coefficients[,2]
    X <- matrix(rep(cofactor, r), nrow = r)
    dat <- dat - beta * X
    # calculate pheno explaine part
    design = model.matrix(~ pheno, dd)
    # fit is introduced design.matrix of pheno
    fit <- lmFit(dat, design);
    fit1 <- lmFit(dat);
    # F-test
    RSS <- rowSums((dat - fitted(fit))^2)
    RSS1 <- rowSums((dat - fitted(fit1))^2)
    df1 <- length(levels(pheno)) + 1 - 1
    df2 = n - length(levels(pheno)) - 1
    F = (RSS1 - RSS) / df1 / (RSS / df2);
    if(df2 > 1e6){
      F.p.value = pchisq(df1 * F, df1, lower.tail = FALSE)
    }
    else{
      F.p.value = pf(F, df1, df2, lower.tail = FALSE)
    }
    tab <- data.frame(intercept=fit$coefficients[,1], beta = fit$coefficients[,2], f = F, pvalue =F.p.value)
  }
  else if(type == "categorical" && is.null(cofactor)){
    cat(sprintf("there is no cofactor"))
    pheno = as.factor(pheno)
    design = model.matrix(~ pheno)
    fit <- lmFit(dat, design)
    fit1 <- lmFit(dat)
    # F-test
    RSS = rowSums((dat - fitted(fit))^2)
    RSS1 = rowSums((dat - fitted(fit1))^2)
    df1 = length(levels(pheno)) - 1;
    df2 = n - length(levels(pheno));
    F = (RSS1 - RSS) / df1 / (RSS / df2);
    if(df2 > 1e6){
      F.p.value = pchisq(df1 * F, df1, lower.tail = FALSE)
    }
    else{
      F.p.value = pf(F, df1, df2, lower.tail = FALSE)
    }
    tab <- data.frame(intercept=fit$coefficients[,1], beta = fit$coefficients[,2], f = F, pvalue =F.p.value)
  }
  else if(type == "continuous" && !is.null(cofactor)){
    cat(sprintf("cofactor is set as continuous variables for continuous pheno ..."))
    pheno = as.numeric(pheno);
    if(is.factor(pheno))
      pheno = .factor2num(pheno)
    else
      pheno = as.numeric(pheno)
    dd = data.frame(pheno = pheno, cofactor = cofactor)
    cofactor.design = model.matrix(~ cofactor, dd)
    fit <- lmFit(dat, cofactor.design)
    # leave cofactor-related part out
    beta = fit$coefficients[,2]
    X <- matrix(rep(cofactor, r), nrow = r)
    dat <- dat - beta * X
    # calculate pheno explaine part
    design = model.matrix(~ pheno, dd)
    # fit is introduced design.matrix of pheno
    fit <- lmFit(dat, design);
    fit1 <- lmFit(dat);
    # F-test
    RSS <- rowSums((dd - fitted(fit))^2)
    RSS1 <- rowSums((dd - fitted(fit1))^2)
    df1 <- length(levels(pheno)) + 1 - 1
    df2 = n - length(levels(pheno)) - 1
    F = (RSS1 - RSS) / df1 / (RSS / df2);
    if(df2 > 1e6){
      F.p.value = pchisq(df1 * F, df1, lower.tail = FALSE)
    }
    else{
      F.p.value = pf(F, df1, df2, lower.tail = FALSE)
    }
    tab <- data.frame(intercept=fit$coefficients[,1], beta = fit$coefficients[,2], f = F, pvalue =F.p.value)
  }
  else if(type == "continuous" && is.null(cofactor)){
    cat(sprintf("there is no cofactor for continuous pheno regression ..."))
    if(is.factor(pheno))
      pheno = .factor2num(pheno)
    else
      pheno = as.numeric(pheno)
    design = model.matrix(~ pheno)
    fit <- lmFit(dat, design)
    fit1 <- lmFit(dat)
    # F-test
    RSS = rowSums((dat - fitted(fit))^2)
    RSS1 = rowSums((dat - fitted(fit1))^2)
    df1 = length(levels(pheno)) - 1;
    df2 = n - length(levels(pheno));
    F = (RSS1 - RSS) / df1 / (RSS / df2);
    if(df2 > 1e6){
      F.p.value = pchisq(df1 * F, df1, lower.tail = FALSE)
    }
    else{
      F.p.value = pf(F, df1, df2, lower.tail = FALSE)
    }
    tab <- data.frame(intercept=fit$coefficients[,1], beta = fit$coefficients[,2], f = F, pvalue =F.p.value)
  }
  # FDR test By the assumption that null label is a uniform distribution
  p0 <- pi0.est(tab$pvalue[!is.na(tab$pvalue)])$p0;
  tab$qvalue <- qvalue.cal(tab$pvalue, p0)
  tab = tab[tab$qvalue <= FDR, ]
  o = order(tab$pvalue)
  tab <- tab[o,]
  tab <- tab[tab$pvalue <= pval, ]
  if(nrow(tab) == 0){
    cat(sprintf("no signifcant DMP within FDR cutoff ... \n"))
  }
  return(tab);
}

## function    : .factor2num 
## description : utility function to convert facotor variabes to x to numeric
## @argument   : x
## @output     : numeric value of x
.factor2num <- function(x) {
  return(as.numeric(levels(x))[x])
}

## function    : DMP_finder
## description : Generate Differential Methylation Positions From BMarray object
## @argument   : object   --- BMarray object
##               FDR      --- false discovery rate
##               pval     --- pvalue
##               attr     --- measurement type c("M", "beta")
##               cofactor --- cofactor attribute's name ~ character
##               type     --- treat BMarray's type classification as categorical factor or continuous factor
## @output     : table    --- significant position table 
DMP_finder <- function(object = NULL, attr = c("M", "beta"), cofactor = c("age", "gender"), type = c("categorical", "continuous"), FDR = 0.05, pval = 0.05) {
  options(warn = -1)
  if(is.null(object) || !"BMarray" %in% class(object))
    stop("BMarray object is necessary to extract attribute's information")
  attr = match.arg(attr)
  cofactor = match.arg(cofactor)
  dat = getBeta(object, type = attr)
  pheno = .groud_id(object, "Type")
  x = .get_cofactor(merge.dataframe)[.groud_id(object, "Pair"),][[cofactor]]
  x = .factor2num(x)
  cat(sprintf("pheno size = %d \n", length(pheno)))
  cat(sprintf("cofactor size = %d \n", length(x)))
  dmp = DMP_engine(dat, pheno = pheno, cofactor = x, type = type, FDR = FDR, pval = pval)
  return(dmp);
}

.unittest_DMP_finder <- function(){
  ## result.unpaired is from .unittest_Assembler; 
  ## load("../data/BMarrayDB_all.rda")
  dmp = DMP_finder(result.unpaired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.05, pval = 0.05)
  dmps = DMP_finder(result.unpaired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.01, pval = 0.05)
  save(dmps, file = "../data/dmp_table.Rda")
}


## DMR search : Based on bumphunter
require(doParallel)
require(bumphunter)

## @Method     : bumphunter
## @argument   : BMarray object and parameters for bumphunterEngine
## @output     : data frame of DMR(differential methylation region)
## description : this method is similar to bumphunter for minif
## @change     : new version bumphunterEngine delete nullmethod; onlu permutation is allowed
setMethod("bumphunter", signature(object = "BMarray"),
          function(object, design, cluster=NULL,
                   coef = 2,
                   cutoff=NULL, pickCutoff=FALSE, pickCutoffQ=0.99, 
                   maxGap = 500,
                   smooth = FALSE,
                   smoothFunction = locfitByCluster,
                   useWeights = FALSE,
                   B=ncol(permutations), permutations=NULL,
                   verbose = TRUE,
                   type = c("beta","M"), ...) {
            type <- match.arg(type)
            bumphunterEngine(getBeta(object, type),
                             design = design,
                             chr = as.factor(seqnames(object)),
                             pos = start(object),
                             cluster = cluster,
                             coef = coef,
                             cutoff=cutoff,
                             pickCutoff=pickCutoff,
                             pickCutoffQ=pickCutoffQ, 
                             maxGap = maxGap,
                             smooth = smooth,
                             smoothFunction = smoothFunction,
                             useWeights = useWeights,
                             B=B,
                             permutations=permutations,
                             verbose = verbose, ...)
          })

## function    : DMR_finder
## description : DMR_finder is a Differential Methylated Region search engine by bumphunter algorithm
## @arguments  : object     --- BMarray
##               coef       --- which col of design matrix will be used to find methylation differential region
##               nullMethod --- To find out significant methylated region, we will create a H_0 hypothesis(null distribution)
##                              Therefore, there are 2 methods to generate null distribution : {"permutation","bootstrap"}
##               smooth     --- whether smooth the Beta or M value of nearest CpG sites [within 1000 nt]
##               useWeights --- smooth should use weights or not
##               mcore      --- multi-thread parallel
##               type       --- measurement to analysis DMR : {"beta", "M"}
##               maxGap     --- max distance of two CpG sites within same cluster
##               B          --- How many times permutation to generation null distribution
## @output     : dmr data.frame
DMR_finder <- function(object, coef = 2,  smooth = TRUE, useWeights = TRUE, mcore = 4, type = c("beta", "M"),
                       maxGap = 500, B = 100) {
  options(warn = -1)
  type = match.arg(type)
  # nullMethod = match.arg(nullMethod)
  pheno <- pData(object)$Type
  designMatrix <- model.matrix(~ pheno)
  registerDoParallel(cores = mcore)
  cat(sprintf("start bumphunting ...\n"))
  # dmrs <- bumphunter(M, design = designMatrix, chr = chr, pos = pos, coef = coef, pickCutoff = TRUE, B=B, smooth=smooth, useWeights = useWeights, maxGap = maxGap)
  dmrs <- bumphunter(object, design = designMatrix, coef = coef, 
                     maxGap = maxGap, pickCutoff = TRUE, B=B, 
                     smooth=smooth, useWeights = useWeights, type = type)
  dmrs
}

.unittest_DMR_finder <- function(){
  dmrs <- DMR_finder(result.unpaired, nullMethod = "permutation", type = "beta", B=100)
}


.unittest_All_finder <- function(){
  cat("loading BMarray Database ...")
  load("../data/BMarrayDB.rda")
  result.paired = tmp$BMarray
  rm(tmp)
  ## doParallel should claim cores to use ahead
  registerDoParallel(cores = 6)
  # 500 nt cluster
  dmrs_gap500  <- DMR_finder(result.paired, mcore = 6, type = "M", B=100)
  dmrs_gap500  <- DMR_annotation(object = result.paired, dmr = dmrs_gap500, type = "M")
  # 1000 nt cluster
  dmrs_gap1000 <- DMR_finder(result.paired, mcore = 6, type = "M", maxGap = 1000, B=100)
  dmrs_gap1000  <- DMR_annotation(object = result.paired, dmr = dmrs_gap1000, type = "M")
  # 2000 nt cluster
  dmrs_gap2000 <- DMR_finder(result.paired, mcore = 6, type = "M", maxGap = 2000, B=100)
  dmrs_gap2000  <- DMR_annotation(object = result.paired, dmr = dmrs_gap2000, type = "M")
  ## save
  save(dmrs_gap500, file = "../data/dmrs_gap500.Rda")
  save(dmrs_gap1000, file = "../data/dmrs_gap1000.Rda")
  save(dmrs_gap2000, file = "../data/dmrs_gap2000.Rda")
  rm(dmrs_gap500);
  rm(dmrs_gap1000);
  rm(dmrs_gap2000);
  
  ## DMP_finder : different FDR
  dmp_fdr.05 = DMP_finder(object = result.paired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.05, pval = 0.05)
  dmp_fdr.01 = DMP_finder(object = result.paired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.01, pval = 0.05)
  ## save
  save(dmp_fdr.05, file = "../data/dmp_fdr5.Rda")
  save(dmp_fdr.01, file = "../data/dmp_fdr1.Rda")
  rm(dmp_fdr.05);
  rm(dmp_fdr.01);
}

.unittest_DMP_all <- function(){
  load("../data/BMarrayDB_all.rda")
  result.unpaired = unpairedDB$BMarray
  dmp_fdr.05 = DMP_finder(object = result.unpaired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.05, pval = 0.05)
  dmp_fdr.01 = DMP_finder(object = result.unpaired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.01, pval = 0.05)
  ## save
  save(dmp_fdr.05, file = "../data/dmp_fdr5_all.Rda")
  save(dmp_fdr.01, file = "../data/dmp_fdr1_all.Rda")
  rm(dmp_fdr.05);
  rm(dmp_fdr.01);
}


.unittest_DMP_pair <- function(){
  cat("loading BMarray Database ...\n")
  load("../data/BMarrayDB.rda")
  result.paired = tmp$BMarray
  rm(tmp)
  ## DMP_finder : different FDR
  dmp_fdr.05 = DMP_finder(object = result.paired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.05, pval = 0.05)
  dmp_fdr.01 = DMP_finder(object = result.paired, attr = "M", cofactor = "age", type = "categorical", FDR = 0.01, pval = 0.05)
  ## save
  save(dmp_fdr.05, file = "../data/dmp_fdr5.Rda")
  save(dmp_fdr.01, file = "../data/dmp_fdr1.Rda")
  rm(dmp_fdr.05);
  rm(dmp_fdr.01);
}






