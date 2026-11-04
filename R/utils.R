# Part of the getBio package.

## get tpm  ##
#' get tpm
#'
#' @export
get_tpm = function(data){
  fpkmToTpm <- function(fpkm)
  {
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  
  tpms <- apply(data,2,fpkmToTpm)
  return(tpms)
}


# get top 'n' genes from an ordered Data #
#' get top 'n' genes from an ordered Data
#'
#' @export
Top = function(data,n){
  n = as.numeric(n)
  upgene = data[1:n,]
  downgene = tail(data,n)
  matrix = rbind(upgene, downgene)
  return(matrix)
}

# get difference analysis from limma algorithn #


## get median number of a column or a row in a data frame ##
#' get median number of a column or a row in a data frame
#'
#' @export
Get_Med = function(data){
  library(stats)
  data = as.numeric(data)
  return(median(data))
}


## return logic of the median ##
#' return logic of the median
#'
#' @export
Get_logic = function(data){
  return(data>Get_Med(data))
}


## check if gene exist ##
#' check if gene exist
#'
#' @export
Chect_gene = function(expdata,gene){
  ifelse(length(grep(row.names(expdata),pattern = ''))==0,paste('this gene is no exist'),
         )
  row.names(expdata)[grep(row.names(expdata),pattern = 'l')]
}


## preprosessing and rank the data##
#' preprosessing and rank the data
#'
#' @export
get_Normalize = function(data,choose_normalized=F){
  data = data[rowSums(data<=0.1)<ncol(data),]
  if(choose_normalized==T){
    data=log2(data+1)
  }
  return(data)
}


## preprosessing and rank the data ##
#' preprosessing and rank the data
#'
#' @export
get_rankMatrix = function(choose_data, group_list, Ctrgroup, Expgroup){
  choose_data = as.data.frame(choose_data)
  all_group=unique(group_list)
  length1=length(which(group_list==Ctrgroup))
  Ctrgroup = which(group_list==Ctrgroup)
  length2=length(which(group_list==Expgroup))
  Expgroup = which(group_list==Expgroup)
  ncols=ncol(choose_data)
  if (length1==1&length2==1){
    choose_data$log2fc=choose_data[,Expgroup]-choose_data[,Ctrgroup]
  }else if(length1==1&length2>1){
    choose_data$log2fc=apply(choose_data[,Expgroup],1,mean)-choose_data[,1]
  }else if(length1>1&length2==1){
    choose_data$log2fc=choose_data[,ncols]-apply(choose_data[,Ctrgroup],1,mean)
  }else{
    choose_data$log2fc=apply(choose_data[,Expgroup],1,mean)-
      apply(choose_data[,Ctrgroup],1,mean)
  }
  choose_data=choose_data[order(choose_data[,'log2fc'],decreasing=TRUE),]
  return(choose_data)
}


## reannotate gene from GEO data only with probe ID ##
#' reannotate gene from GEO data only with probe ID
#'
#' @export
get_reannotated = function(gpl,id,exprset,method = 'mean'){
  library(idmap3)
  probe2gene=idmap3::get_pipe_IDs(gpl)
  
  library(dplyr)
  exprset = cbind(id,exprset)
  colnames(probe2gene) = c(colnames(exprset)[1],"symbol")
  exprset <- merge(probe2gene,exprset,by='id')
  exprset <- exprset%>%dplyr::arrange(symbol)
  ifelse(method == 'mean',
         # use mean as certain gene expression level
         exprset2 <- aggregate(x = exprset[,3:ncol(exprset)],
                               by = list(exprset$symbol),
                               FUN = mean),
         # use medianas certain gene expression level
         exprset2 <- aggregate(x = exprset[,3:ncol(exprset)],
                               by = list(exprset$symbol),
                               FUN = median)
  )
  return(exprset2)
}


## get the numic format matrix of expression data frame ##
#' get the numic format matrix of expression data frame
#'
#' @export
get_numMatrix = function(data){
  data=as.matrix(data)
  dimnames=list(rownames(data),colnames(data))
  data=matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
  data=avereps(data)
  return(data)
}


## figure out if one gene (row name) exist in as data frame
#' figure out if one gene (row name) exist in as data frame
#'
#' @export
get_genes = function(data,pattern){
  return(row.names(data)[grep(row.names(data),pattern = pattern)])
}


## find the correlation genes of certain given gene ##
#' find the correlation genes of certain given gene
#'
#' @export
get_correlated_genes <- function(seurat_object = NULL,expression_matrix = NULL, gene_of_interest, num_genes = 50) {
  if (is.null(seurat_object)) {
    gene_expression <- expression_matrix[gene_of_interest, ]
    correlation_values <- apply(expression_matrix, 1, function(x) cor(x, gene_expression))
    sorted_genes <- sort(correlation_values, decreasing = TRUE)
    top_genes <- sorted_genes[1:num_genes]
    return(top_genes)
  }else{
    message("find correlation genes of any given gene, scRNAseq only")
    gene_expression <- GetAssayData(seurat_object, slot = "data")[gene_of_interest,]
    correlation_values <- apply(GetAssayData(seurat_object, slot = "data"), 1, function(x) cor(x, gene_expression))
    sorted_genes <- sort(correlation_values, decreasing = TRUE)
    top_genes <- sorted_genes[1:num_genes]
    return(top_genes)
  }
  
}
