## Annotation script is for DMR region's annotation and mapping our DMR on hg19 reference
## genomics
## In addition, we introduce DMR class via Annotation.R

require(BiocGenerics)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## S4 object   :  DMR class constructor
setClass("DMRs", contains = "GRanges",
         representation(annotation = "character", dmrs = "data.frame"))

## function    :  DMR_annotation
## description :  DMR constructor
## @argument   :  object  --- BMarray object : information retrieval handel
##                dmrs    --- dmr(differential methylation region) table(data.frame)
##                type    --- c("beta", "M")
## @output     :  DMR object
DMR_annotation <- function(object = NULL, dmr = NULL, type = c("beta", "M")) {
  options(warn = -1)
  if(is.null(object) || is.null(dmr))
    stop("BMarray object and dmr data frame are needed ...")
  stopifnot("BMarray" %in% class(object))
  # methylation measurement type
  type = match.arg(type)
  # human genome reference's annotation From IlluminaHumanMethylation450kanno.ilmn12.hg19
  annotation = getAnnotation(object)
  mat = getBeta(object, type)
  # get DMR$table : DMR is result of bumphunter
  dmrs <- dmr$table
  # pdmr -- parameters of dmr
  pdmr <- as.data.frame(t(apply(dmrs, 1, function(dat){
    chr <- dat[1]; 
    start <- as.numeric(dat[2]); end <- as.numeric(dat[3]);
    beta0 <- as.numeric(dat[4]);
    cgs <- paste(rownames(mat)[as.numeric(dat[7]) : as.numeric(dat[8])], collapse = ";")
    Islands_Name <- paste(unique(annotation@data$Islands.UCSC[rownames(mat)[as.numeric(dat[7]) : as.numeric(dat[8])], ][,"Islands_Name"]), 
                          collapse = ";")
    Relation_to_Island <- paste(unique(annotation@data$Islands.UCSC[rownames(mat)[as.numeric(dat[7]) : as.numeric(dat[8])], ][,"Relation_to_Island"]), 
                                collapse = ";")
    pvalue <- as.numeric(dat[13])
    return(c(chr, start, end, beta0, cgs, Islands_Name, Relation_to_Island, pvalue))
  })), stringsAsFactors=FALSE)
  # name of pdmr
  names(pdmr) <- c("chr", "start", "end", "beta0", "cgs", "Islands_Name", "Relation_to_Island", "pvalue")
  DMRs <- GRanges(seqnames = pdmr$chr, ranges = IRanges(start = as.numeric(pdmr$start), end = as.numeric(pdmr$end)), 
                  beta0 = as.numeric(pdmr$beta0), cgs = Rle(pdmr$cgs), Islands_Name = Rle(pdmr$Islands_Name), 
                  Relation_to_Island = Rle(pdmr$Relation_to_Island), pvalue = as.numeric(pdmr$pvalue));
  DMRs <- as(DMRs, "DMRs")
  DMRs@annotation <- unlist(list(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))
  DMRs@dmrs <- as.data.frame(dmrs);
  DMRs 
}


.unittest_annotation <- function(){
  aDMRs <- DMR_annotation(object = result.unpaired, dmr = dmrs, type = "M")
}

