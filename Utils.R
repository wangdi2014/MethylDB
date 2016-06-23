## Utils.R is a script which collect ultility function which perhaps used for Data analysis
## warning : methyl450kpatient is for 450K illumina data information
##           merge.dataframe of "BRCA_metadb.rda" is patients' barcode, gender and age information
require(MASS)

## function    : .group_id 
## description : return BMarray object's patient group information; data type is `tumor` or `normal`
## @argument   : object  -- should be validate as BMarray object
## @output     : type vector
.groud_id <- function(object, type) {
  stopifnot("BMarray" %in% class(object))
  return(pData(object)[,type])
  #return(pData(object)[,"Type"]);
}

## function    : .get_cofactor
## description : return a age related dataframe of all the other cofactor
## @arguments  : merge.dataframe (specific data)
## @output     : dataframe of merge.dataframe but modify the original TCGA-barcode; just keep first 3 chars
##               e.g ==> row TCGA-E2-A1BD-01A-11D-A12R-05 's name is TCGA-E2-A1BD
.get_cofactor <- function(df){
  names <- as.character(df$tcga_idx);
  names <- unlist(lapply(str_split(names, "-"), function(x){return(paste(x[1:3], collapse="-"))}));
  df$tcga_idx <- names;
  df <- unique(df)
  rownames(df) <- df$tcga_idx
  return(as.data.frame(df));
}


## function    : .write_matrix
## description : write a edge based matrix into 
## @arguments  : edge_mat   --- edge table
##               filename   --- output file name
## @output     : file
.write_matrix <- function(edge_mat = NULL, filename = NULL) {
  options(warn = -1)
  if(is.null(edge_mat) || is.null(filename))
    stop("edge matrix and filename are needed ...")
  write.table(edge_mat, file = filename, col.names = FALSE, row.names = FALSE)
}
