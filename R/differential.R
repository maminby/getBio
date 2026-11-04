# Part of the getBio package.

#' Get_design
#'
#' @export
Get_design = function(tpms, group_list, pairinfo = NULL){
  #pairinfo = factor(rep(1:nrow(group1),number(groups)))
  library(clusterProfiler)
  library(DOSE)
  if (is.null(pairinfo)) {
    group_list <- factor(group_list,levels = unique(group_list),ordered = F)
    
    # difference analysis
    design <- model.matrix(~0+group_list)
    colnames(design) = levels(group_list)
    rownames(design) = colnames(tpms)
    return(design)
  }else{
    group_list <- factor(group_list,levels = unique(group_list),ordered = F)
    design <- model.matrix(~group_list+pairinfo)
    rownames(design) = colnames(tpms)
    return(design)
  }
  
}


#' get_diffecence
#'
#' @export
get_diffecence = function(tpms,design,contrast.matrix,logFC_t = 1,pair = F,Pvalue = 0.05){
  if(isTRUE(pair)){
    ## step1 ##
    fit <- lmFit(tpms,design)
    
    ## step2 ##
    fit2 <- eBayes(fit)
    
    ## step3 ##
    tempOutput = topTable(fit2,adjust='fdr',coef=2,number=Inf ,p.value=0.05) 
    nrDEG = na.omit(tempOutput)
    nrDEG$g = ifelse(nrDEG$P.Value>Pvalue,'stable',
                     ifelse( nrDEG$logFC > logFC_t,'UP',
                             ifelse( nrDEG$logFC < -logFC_t,'DOWN','stable')))
  }else{
    ## step1 ##
    fit <- lmFit(tpms,design)
    
    ## step2 ##
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    ## step3 ##
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput)
    
    
    nrDEG$g = ifelse(nrDEG$P.Value>Pvalue,'stable',
                     ifelse( nrDEG$logFC > logFC_t,'UP',
                             ifelse( nrDEG$logFC < -logFC_t,'DOWN','stable')))
  }
  nrDEG$symbol = rownames(nrDEG)
  return(nrDEG)
}


#' run_limmatrend
#'
#' @export
run_limmatrend <- function(L) {
  library(limma)
  library(edgeR)
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(.get_assay_data(L, "counts"), group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~L$condt)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(y, design = design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
  })
  
  hist(tt$P.Value, 50)
  hist(tt$adj.P.Val, 50)
  limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
  plotMD(fit)
  
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$P.Value,
                       padj = tt$adj.P.Val,
                       row.names = rownames(tt)))
}


#' run_voomlimma
#'
#' @export
run_voomlimma <- function(L) {
  library(limma)
  library(edgeR)
  message("voomlimma")
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(L$count, group = L$condt)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~L$condt)
    vm <- voom(dge, design = design, plot = TRUE)
    fit <- lmFit(vm, design = design)
    fit <- eBayes(fit)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
  })
  
  hist(tt$P.Value, 50)
  hist(tt$adj.P.Val, 50)
  limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
  plotMD(fit)
  
  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$P.Value,
                       padj = tt$adj.P.Val,
                       row.names = rownames(tt)))
}
