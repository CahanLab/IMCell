########################################################################
#
#	Useful functions
#
########################################################################



#' loads an R object when you don't know the name
#'
#' loads an R object when you don't know the name
#' @param fname file
#'
#' @return variable
#'
#' @export
utils_loadObject<-function
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}


#' Loads .obs and .var data from loom file
#'
#' @param path path to loom file
#' @param col_attrs column attributes to include in sample table. By default will include all in loom file.
#' @param obs_names where in loom file observation names are stored
#' @param var_names where in loom file variable names are stored
#'
#' @return list includeing expression matrix and sample table
#' 
#' @export
#'
loadDataFromLoom<-function(path,
                    col_attrs=NULL,
                    obs_names='obs_names',
                    var_names='var_names'
  ){

  lfile <- connect(filename = path, skip.validate=TRUE)

  # set obs_names and var_names
  if ("cell_names" %in% lfile[['col_attrs']]$names){
    obs_names<-"cell_names"
  }
  if ("gene_names" %in% lfile[['row_attrs']]$names){
    var_names<-"gene_names"
  }

  # get column attributes to extract
  if (is.null(col_attrs)){
    col_attrs<-lfile[['col_attrs']]$names
  }
  if (sum(!(col_attrs %in% lfile[['col_attrs']]$names))>0){
    stop("elements in col_attrs not in metadata.")
  }

  sampTab<-tryCatch({lfile$get.attribute.df(attributes = col_attrs, col.names=obs_names)},
  			error=function(e){
  				lfile$get.attribute.df(attribute.names = col_attrs, col.names=obs_names)
  			})

  geneNames<-lfile[["row_attrs"]][[var_names]][]
  cellNames<-lfile[["col_attrs"]][[obs_names]][]
  expMat<- t(lfile[["matrix"]][,])
  rownames(expMat)<-geneNames
  colnames(expMat)<-cellNames

  list(expDat = expMat, sampTab = sampTab)
}



#' Downsamples and log-transform counts data
#'
#' @param counts raw counts expression matrix
#' @param xFact what to dowmsample and normalize to
#' @param observations_by_features row and column orientation of the counts matrix
#'
#' @return normalized and transformed expression matrix
#' 
#' @export
#'
downsample_and_transform<-function(counts,xFact=1e5,observations_by_features=FALSE){
	if (!observations_by_features){
		counts<-t(counts)
	}
	rowsum<-rowSums(counts)
	rowsum[rowsum==0]<-1
	expX<-sweep(counts,1,rowsum,FUN='/')
	expX<-expX*xFact
	expX<-log2(expX+1)

	if (!observations_by_features){
		expX<-t(expX)
	}
	
	expX
}








