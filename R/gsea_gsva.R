# Part of the getBio package.

#' get_gsea
#'
#' @export
get_gsea = function(data,datatype = 'nrDEG',geneset,genesetNum,filename,height = 5,width = 7){
  #require('RColorBrewer')
  #colors <- colorRampPalette( rev(brewer.pal(9, "RdYlBu")))(255)
  if (datatype == 'nrDEG') {
    # this function only used for nrDEG
    genelist = data[,c(8,3)]
    de1 <- as.numeric(genelist[,2])
    names(de1) = genelist[,1]
    de1 = sort(de1,decreasing = T)
    
    y <- GSEA(de1, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff = 1)
    p <- gseaplot2(y, 1:genesetNum,base_size=10,subplots = 1:4,title = filename,pvalue_table = T)
    graph2ppt(x=p, file = filename, height = height,width = width, aspectr=4)
  }
  else{
    # this function only used for self-defined difference expression matrix
    genelist = data
    de1 <- as.numeric(genelist[,2])
    names(de1) = genelist[,1]
    de1 = sort(de1,decreasing = T)
    
    y <- GSEA(de1, TERM2GENE=geneset, verbose=FALSE,pvalueCutoff = 1)
    p <- gseaplot2(y, 1:genesetNum,base_size=10,subplots = 1:4,title = filename,pvalue_table = T)
    graph2ppt(x=p, file = filename, height = height,width = width, aspectr=4)
  }
}


#' get_gsva
#'
#' @export
get_gsva = function(data, gene_set, gene_setname, species = 'm',save =F){
  if (!suppressWarnings(suppressMessages(require("GSVA", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    install.packages(p, repos=site)
    require("GSVA",character.only=T)
  }
  if (species=='m'|species=='mouse') {
    gene_set = str_to_title(gene_set)
  }else {
    gene_set = str_to_upper(gene_set)
  }
  
  list = split(as.matrix(gene_set)[,1], rep(gene_setname,length(gene_set)))
  data = get_numMatrix(data)
  
  ssgseaScore=gsva(data, list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  gsva_matrix1<- t(scale(t(ssgseaScore)))
  return(gsva_matrix1)
  if (save == T) {
    write.table(gsva_matrix1,file=paste0(gene_setname,"ssgseaOut_totalgene_scale.txt"),sep="\t",quote=F)
  }
}


#' get_immune_gsva
#'
#' @export
get_immune_gsva = function(data, group_list, gene,save = F){
  library(pheatmap)
  if (!suppressWarnings(suppressMessages(require("GSVA", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    install.packages(p, repos=site)
    require("GSVA",character.only=T)
  }
  gene_set<-read.csv("../mmc3.csv")[, 1:2]
  gene_set$Metagene = str_to_title(gene_set$Metagene)
  
  list<- split(as.matrix(gene_set)[,1], gene_set[,2])
  data = get_numMatrix(data)
  ssgseaScore=gsva(data, list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
  gsva_matrix<-ssgseaScore
  
  
  #缩放
  gsva_matrix1<- t(scale(t(gsva_matrix)))
  write.table(gsva_matrix1,file="ssgseaOut_totalgene_scale.txt",sep="\t",quote=F)
  
  anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
  pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
  anti<- gsub('^ ','',rownames(gsva_matrix1))%in%anti_tumor #gsub是替代 %in%是询问前者是否为后者
  pro<- gsub('^ ','',rownames(gsva_matrix1))%in%pro_tumor
  non <- !(anti|pro)    #既不是anti或pro的赋值为non
  
  #将anti pro non 分开合并
  gsva_matrix1<- rbind(gsva_matrix1[anti,],gsva_matrix1[pro,],gsva_matrix1[non,])
  #定义ssGSEA score矫正函数
  normalization<-function(x){
    return((x-min(x))/(max(x)-min(x)))}
  nor_gsva_matrix1 <- normalization(gsva_matrix1)
  write.table(nor_gsva_matrix1,file="ssgseaOut_mmc3.txt",sep="\t",quote=F)
  annotation_col = data.frame('Group'=factor(group_list,levels = unique(group_list)))
  rownames(annotation_col)<-colnames(nor_gsva_matrix1)
  p = pheatmap(nor_gsva_matrix1[,order(data[gene,],decreasing = F)],
               show_colnames = F,
               cluster_rows = F,cluster_cols = F,
               annotation_col = annotation_col)
  
  print(p)
  
  if (save == T) {
    pdf(paste(gene,'immune_GSVA.pdf',sep = '_'),width = 10,height = 10)
    print(p)
    dev.off()
  }
  
}


#' get_broadinstituteORG
#'
#' @export
get_broadinstituteORG = function(path){
  message('path should be upper folder of "cluster","metadata"and "expression"')
  #cluster.path = list.files(paste0(path,'/cluster'),full.names = T)
  mtx.path = list.files(list.dirs(paste0(path,'/expression'))[2],pattern = 'mtx',full.names = T)
  cells.path = list.files(list.dirs(paste0(path,'/expression'))[2],pattern = 'Barcodes',full.names = T)
  features.path = list.files(list.dirs(paste0(path,'/expression'))[2],pattern = 'Genes',full.names = T)
  metadata.path = list.files(paste0(path,'/metadata'),full.names = T)
  
  ## read in expression ##
  expression = ReadMtx(mtx = mtx.path,cells = cells.path,features = features.path)
  expression = Seurat::CreateSeuratObject(expression)
  expression = get_newReduction(expression)
  ## add metadata ##
  test = read.csv(cluster.path,header = T,row.names = 1)
  expression@meta.data = read.csv(metadata.path,header = T,row.names = 1)[-1,]
  return(expression)
}


#' GeneSet S4 class
#'
#' @slot gene_set_name,gene_set_description,gene_names character
#' @exportClass GeneSet
#' @importFrom methods new
##  YOU SHALL ANDERSTAND WHAT THIS EXACTLY DO!  ##
# class create
methods::setClass("GeneSet",
         methods::representation(
           gene_set_name = "character",
           gene_set_description = "character",
           gene_names = "character"
         )
)


#' GeneSet
#'
#' @export
GeneSet <- function(gene_set_name, gene_set_description, gene_names) {
  new("GeneSet", gene_set_name = gene_set_name, gene_set_description = gene_set_description, gene_names = gene_names)
}


#' createGeneSetCollection
#'
#' @export
createGeneSetCollection <- function(gmt_file) {
  gmt_lines <- readLines(gmt_file)
  gene_sets <- list()
  for (line in gmt_lines) {
    elements <- strsplit(line, "\t")[[1]]
    gene_set_name <- elements[1]
    gene_set_description <- elements[2]
    gene_names <- elements[3:length(elements)]
    gene_set <- GeneSet(gene_set_name, gene_set_description, gene_names)
    gene_sets[[gene_set_name]] <- gene_set
  }
  return(gene_sets)
}
