# Part of the getBio package.

#### get the peso cell from certain cluster ####
#' get the peso cell from certain cluster
#'
#' @export
GatherData.Seurat <- function(object,
                              assay,
                              slot_use,
                              ...) {
  
  assay <- assay %||% "RNA"
  slot_use <- slot_use %||% "data"
  obj_data <- .get_assay_data(object, what = slot_use, assay = assay) %>%
    as.matrix()
  return(obj_data)
}


#' PseudoCell
#'
#' @export
PseudoCell  <- function(object,
                        assay_use = NULL,
                        slot_use = NULL,
                        cluster_use =NULL,
                        pseudocell.size  =NULL){
  message("tips:
          object : count of seurat object
          Cluster_use : one col in metadata
          pseudocell.size : how many cell will be pseudo")
  
  Inter<- GatherData.Seurat(object = object,
                     assay = assay_use,
                     slot_use = slot_use) 
  Inter[Inter<0]=0
  idd<-object@meta.data
  Inter.id<-cbind(rownames(idd),as.vector(idd[,cluster_use]))
  
  rownames(Inter.id)<-rownames(idd)
  colnames(Inter.id)<-c("CellID","Celltype")
  
  Inter.id<-as.data.frame(Inter.id)
  Inter1<-Inter[,Inter.id$CellID]
  Inter<-as.matrix(Inter1)
  pseudocell.size = pseudocell.size ## 10 test
  new_ids_list = list()
  Inter.id$Celltype <- as.factor(Inter.id$Celltype)
  for (i in 1:length(levels(Inter.id$Celltype))) {
    cluster_id = levels(Inter.id$Celltype)[i]
    cluster_cells <- rownames(Inter.id[Inter.id$Celltype == cluster_id,])
    cluster_size <- length(cluster_cells)       
    pseudo_ids <- floor(seq_along(cluster_cells)/pseudocell.size)
    pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
    names(pseudo_ids) <- sample(cluster_cells)  
    new_ids_list[[i]] <- pseudo_ids     
  }
  
  new_ids <- unlist(new_ids_list)
  new_ids <- as.data.frame(new_ids)
  new_ids_length <- table(new_ids)
  
  new_colnames <- rownames(new_ids)  ###add
  all.data<-Inter[,as.character(new_colnames)] ###add
  all.data <- t(all.data)###add
  
  new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                      list(name=new_ids[,1]),FUN=mean)
  rownames(new.data)<-new.data$name
  new.data<-new.data[,-1]
  
  new_ids_length<-as.matrix(new_ids_length)##
  short<-which(new_ids_length< pseudocell.size -1 )##
  new_good_ids<-as.matrix(new_ids_length[-short,])##
  result<-t(new.data)[,rownames(new_good_ids)]
  rownames(result)<-rownames(Inter)
  
  newobject <- CreateSeuratObject(result)
  newobject@misc$idtrans <- new_ids
  newobject@commands$PseudoCell <- LogSeuratCommand(newobject, return.command = TRUE)
  return(newobject)
}


#get the cell chat circle plot #
#' get the cell chat circle plot
#'
#' @export
get_Cellchat = function(object,ident.use = "orig.ident"){
  library(CellChat)
  library(patchwork)
  cellchat <- createCellChat(object = object, meta = object@meta.data,group.by = "orig.ident")
  cellchat <- addMeta(cellchat, meta = object@meta.data)
  cellchat <- setIdent(cellchat, ident.use = ident.use) # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mar = rep(1,4),mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  return(cellchat)
}


#get the top gene heatmap of different cluster #
#' get the top gene heatmap of different cluster
#'
#' @export
get_Topmarkers=function(seurat_object,filename){
  library(dplyr)
  library(stringr)
  all.markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
  top5 <- all.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  top10 <- all.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  top20=as.data.frame(top20)
  write.csv(top20,paste0(filename,"_clustertop20_gene.csv"),row.names = F)
  pdf(paste0("heatmap_",filename,"_top5.pdf"),width =15, height = 10.3)
  f1=DoHeatmap(object = seurat_object, features = top5$gene,label = TRUE,size=3)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f1)
  dev.off()
  
  png(paste0("heatmap_",filename,"_top5.png"),width = 1500, height = 980)
  print(f1)
  dev.off()
  
  pdf(paste0("Dotplot_",filename,"_nolegend.pdf"),width =10, height = 5.18)
  f2=DotPlot(object = seurat_object, features = rev(x = unique(top5$gene)), cols = c("blue","red"), 
          dot.scale = 5) + RotatedAxis()+NoLegend()
  print(f2)
  dev.off()
  
}


#' modify_vlnplot
#'
#' @export
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size,cols = color_used,group.by = 'seurat_clusters') +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}


## main function
#' main function
#'
#' @export
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45,vjust = 1,hjust = 1), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


## get new reduction of a seurat object, can be based on certain genes                    
#' get new reduction of a seurat object, can be based on certain genes
#'
#' @export
get_newReduction = function(seurat.object, choose_gene=NULL){
  library(purrr)
  message('choose_gene is the reference genes used for reduction')
  seurat.object = NormalizeData(seurat.object,verbose = T)%>%
    FindVariableFeatures(.,selection.method = "vst")%>%
    ScaleData()
  if (is.null(choose_gene)) {
    seurat.object = RunPCA(seurat.object, verbose = FALSE)
  }
  else{
    seurat.object = RunPCA(seurat.object, verbose = FALSE, features = choose_gene)
  }
  ndims = ifelse(length(Stdev(seurat.object, reduction = 'pca'))>30,
                 30,length(Stdev(seurat.object, reduction = 'pca')))
  seurat.object = FindNeighbors(object = seurat.object, dims = 1:ndims)
  seurat.object = FindClusters(object = seurat.object, resolution = 0.3)%>%
    RunUMAP(., dims = 1:ndims)
  return(seurat.object)
}


## for cluster not split well
#' for cluster not split well
#'
#' @export
get_localReduction = function(Seurat.object, choose_cluster, choose_gene = NULL){
  message('local reducution cluster will store in "new_cluster" slot')
  id.cluster = as.data.frame(Idents(Seurat.object))
  if (length(choose_cluster)>1) {
    for (i in choose_cluster) {
      sub.seurat = subset(Seurat.object,idents = i)
      sub.seurat = get_newReduction(sub.seurat,choose_gene = choose_gene)
      id.cluster[,1] = as.character(id.cluster[,1])
      id.cluster[colnames(sub.seurat),] = paste(i,as.character(Idents(sub.seurat)),sep = '_')
      Seurat.object[['new_cluster']] = factor(id.cluster[,1],levels = unique(id.cluster[,1]))
    }
  }else{
    sub.seurat = subset(Seurat.object,idents = choose_cluster)
    sub.seurat = get_newReduction(sub.seurat,choose_gene = choose_gene)
    id.cluster[,1] = as.character(id.cluster[,1])
    id.cluster[colnames(sub.seurat),] = paste(choose_cluster,as.character(Idents(sub.seurat)),sep = '_')
    Seurat.object[['new_cluster']] = factor(id.cluster[,1],levels = unique(id.cluster[,1]))
  }
  
  return(Seurat.object)
}


## do heatmap of top5 difference genes
#' do heatmap of top5 difference genes
#'
#' @export
myfindmarkers=function(seurat_object,filename){
  library(dplyr)
  library(stringr)
  all.markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
  top5 <- all.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  top10 <- all.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  top5=as.data.frame(top5)
  write.csv(top20,paste0(filename,"_clustertop20_gene.csv"),row.names = F)
  pdf(paste0("heatmap_",filename,"_top5.pdf"),width =15, height = 10.3)
  f1=DoHeatmap(object = seurat_object, features = top5$gene,label = TRUE,size=3,group.colors = color_used)+
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
  print(f1)
  dev.off()
  
  png(paste0("heatmap_",filename,"_top5.png"),width = 1500, height = 980)
  print(f1)
  dev.off()
  
  # pdf(paste0("Dotplot_",filename,"_nolegend.pdf"),width =20, height = 5.18)
  pdf(paste0("Dotplot_",filename,"_nolegend.pdf"),width =6, height = 5.18)
  f2=DotPlot(object = seurat_object, features = rev(x = unique(top5$gene)), cols = c("blue","red")) + 
    RotatedAxis()+NoLegend()
  print(f2)
  dev.off()
}


## replace the idents of origin seurat with new seurat idents 
#' replace the idents of origin seurat with new seurat idents
#'
#' @export
get_Identreplace = function(origin_obj,replace_obj)
{
  message('replace_obj should be included by orgin_obj')
  contain.cell = intersect(colnames(origin_obj), colnames(replace_obj))
  replace_obj = replace_obj[,contain.cell]
  origin_meta = data.frame('cell' = row.names(origin_obj@meta.data),'CellType' = as.character(Idents(origin_obj)))
  replace_meta = data.frame('cell' = row.names(replace_obj@meta.data),'CellType' = as.character(Idents(replace_obj)))
  
  merged_df <- merge(origin_meta, replace_meta, by = "cell", all.x = TRUE)
  merged_df$CellType <- ifelse(is.na(merged_df$CellType.y), as.character(merged_df$CellType.x), merged_df$CellType.y)
  row.names(merged_df) = merged_df$cell
  merged_df = merged_df[colnames(origin_obj),]
  merged_df = merged_df[,c(1,4)]
  Idents(origin_obj) = merged_df$CellType
  return(origin_obj)
}


#' Get_cytometryiod
#'
#' @export
Get_cytometryiod = function(mat.object=NULL,seurat_object=NULL,log=F, genes=c('x','y'),filename=NULL,type=NULL){
  library(ggExtra)
  if (is.null(mat.object)) {
    
    if (is.null(seurat_object)) {
      
      stop("No correct expression matrix supplied!")
      
    }else{
      
      x= genes[1]
      y = genes[2]
      list1 = .get_assay_data(seurat_object, "data")[genes,]
      
      list1 = t(as.matrix(list1))
      list1 = as.data.frame(list1)
      
      xlab = colnames(list1)[1]
      ylab = colnames(list1)[2]
      colnames(list1)[1]= 'x'
      colnames(list1)[2]= 'y'
      
      if (!is.null(type)) {
        p <- ggplot(data = list1,aes(x = x , y = y)) +
          geom_point(color='#9e4aa8') +
          theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
          theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
          theme_classic()+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
          xlab(paste(xlab)) + ylab(paste(ylab)) + geom_xsidedensity(fill = 'lightskyblue') + geom_ysidedensity(fill = 'orangered')
      }else{
        p <- ggplot(data = list1,aes(x = x , y = y)) +
          geom_point(color='#9e4aa8') +
          theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
          theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
          theme_classic()+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
          xlab(paste(xlab)) + ylab(paste(ylab)) + geom_xsidedensity(fill = 'lightskyblue') + geom_ysidedensity(fill = 'orangered')
      }
      
      pdf(paste0("scatterplot_",filename,".pdf"),width =8, height = 8)
      print(p)
      dev.off()
      
      p <- p + stat_cor(data=list1, method = "spearman")
      pdf(paste0("scatterplot_spearman",filename,".pdf"),width =8, height = 8)
      print(p)
      dev.off()
      
    }
    
  }else{
    
    if (log==F) {mat.object = log2(mat.object+1)}
    
    x= genes[1]
    y = genes[2]
    
    list1 = mat.object[row.names(mat.object)%in%genes,]
    
    list1 = t(as.matrix(list1))
    list1 = as.data.frame(list1)
    
    xlab = colnames(list1)[1]
    ylab = colnames(list1)[2]
    colnames(list1)[1]= 'x'
    colnames(list1)[2]= 'y'
    
    
    p <- ggplot(data = list1,aes(x = x , y = y)) +
      geom_point(color='#9e4aa8') +
      theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
      theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
      theme_classic()+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
      xlab(paste(xlab)) + ylab(paste(ylab)) + geom_xsidedensity(fill = 'lightskyblue') + geom_ysidedensity(fill = 'orangered')
    
    pdf(paste0("scatterplot_",filename,".pdf"),width =4, height = 4)
    print(p)
    dev.off()
    
    p <- p + stat_cor(data=list1, method = "spearman")
    pdf(paste0("scatterplot_spearman",filename,".pdf"),width =4, height = 4)
    print(p)
    dev.off()
  }
  
}


## pseudo time analysis ##
#' pseudo time analysis
#'
#' @export
get_pseudotime = function(pesudo){
  library(monocle3)
  data = pesudo
  expression_matrix = .get_assay_data(pesudo, "counts")
  
  cell_metadata <- pesudo@meta.data
  
  gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), stringsAsFactors = F)
  rownames(gene_annotation) <- rownames(expression_matrix)
  
  pesudo = new_cell_data_set(expression_data = expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)
  pesudo = preprocess_cds(pesudo,num_dim = 50)
  pesudo = reduce_dimension(pesudo,preprocess_method = 'PCA')
  
  
  #switch umap embedding graph to original umap
  pesudo = cluster_cells(pesudo,reduction_method = 'UMAP')
  # pesudo@int_colData$reducedDims$UMAP = Embeddings(Treg.bind.TCR,reduction = 'umap')[colnames(pesudo),]
  pesudo@int_colData$reducedDims$UMAP = Embeddings(data,reduction = 'umap')[colnames(pesudo),]
  
  pesudo = learn_graph(pesudo)
  plot_cells(pesudo,group_label_size = 1,graph_label_size = 3)
  
  pseudo <- order_cells(pesudo,reduction_method = 'UMAP')
  return(pseudo)
}


## get harmony for debatch##
#' get harmony for debatch
#'
#' @export
get_harmony = function(seurat){
  library(harmony)
  metadata = seurat@meta.data
  seurat = ScaleData(seurat)
  seurat = seurat %>% RunHarmony("sample", plot_convergence = TRUE)
  seurat = seurat %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
  seurat@meta.data = metadata
  return(seurat)
}


## get seurat go and kegg analysis ##
#' get seurat go and kegg analysis
#'
#' @export
get_seuratGO = function(seurat,ident1 = NULL,
                        ident2 = NULL,species = 'm',
                        filename = NULL){
  if (is.null(slot)|is.null(ident1)|is.null(ident2)) {
    stop("parameters 'slot' or 'ident1/2' is null")
  }
  if (length(unique(Idents(seurat)))>2) {
    seurat = subset(seurat,idents = c(ident1,ident2))
  }
  
  all.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
  upgene = all.markers[all.markers$cluster==ident1,c('gene','p_val_adj')]
  upgene = upgene[order(upgene$p_val_adj,decreasing = F),1]%>%as.character()
  downgene = all.markers[all.markers$cluster==ident2,c('gene','p_val_adj')]
  downgene= downgene[order(downgene$p_val_adj,decreasing = F),1]%>%as.character()
  if (species == 'human'|species=='h') {
    upgene = str_to_upper(upgene)
    downgene = str_to_upper(downgene)
  }else{
    upgene = str_to_title(upgene)
    downgene = str_to_title(downgene)
  }
  return(list('upgene' = upgene,'downgene' = downgene))
}


## read in file after 10X HTO alignment (antibody captured included) ##
#' read in file after 10X HTO alignment (antibody captured included)
#'
#' @export
get_HTOassay = function(object){
  data = CreateSeuratObject(object$`Gene Expression`)
  data[["HTO"]] <- CreateAssayObject(counts = object$`Antibody Capture`)
  data <- NormalizeData(data, assay = "HTO", normalization.method = "CLR")
  return(data)
}


## get the percentage of two slot from Seurat ##
#' get the percentage of two slot from Seurat
#'
#' @export
get_percentage = function(seurat,slot1,slot2,filename = NULL){
  message('slot2 is x axis')
  library(reshape2)
  percent = as.data.frame(melt(table(seurat@meta.data[,slot1],seurat@meta.data[,slot2])))
  ggplot(percent, aes(x = as.factor(Var1), y = value,fill = Var2)) +
    geom_bar(stat = "identity", position = "fill",color = "black", lwd = 0.25) + 
    scale_fill_manual(name = "Clonotype Group",values = c(color_crc)) +
    xlab(label = slot2) + 
    scale_y_continuous(expand = c(0,0))+theme_classic()+
    ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
  ggplot2::ggsave(paste0(filename,'_percent.pdf'),width = 6,height = 4)
}

## create a Seurat object according to data download from single cell broad institute website #


## another myfind marker created by Long J ##
#' another myfind marker created by Long J
#'
#' @export
myfindmarkers2 = function(seurat_object,filename,logfc.threshold=0.25,min.pct=0.25,species="human",all.markers=NULL,colors=colors,
                          fill_colors=colorRampPalette(c("lightskyblue", "white", "red"))(100)){
  if(is.null(all.markers)){
    all.markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, min.pct = min.pct, 
                                  logfc.threshold = logfc.threshold)
  }
  library(dplyr)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, -p_val_adj)
  top50 <- all.markers %>% group_by(cluster) %>% top_n(50, -p_val_adj)
  write.csv(top20,paste0(filename,"_clustertop20_gene.csv"),row.names = F)
  write.csv(top50,paste0(filename,"_clustertop50_gene.csv"),row.names = F)
  all.markers$p_val_adj[which(all.markers$p_val_adj<1e-300)]=1e-300
  all.markers$p_val_adj=-all.markers$p_val_adj
  all.markers$p_val_adj=all.markers$p_val_adj+length(all.markers$p_val_adj):1
  
  if(species=="huamn"){
    common_markers=c("CD3D","CD3E","CD8A","CD8B","CD4")
  }else{
    common_markers=c("Cd3d","Cd3e","Cd8a","Cd8b1","Cd4")
  }
  
  if(length(grep(pattern = "TotalSeqC",rownames(seurat_object)))>0){
    antibody_name=rownames(seurat_object)[grep(pattern = "TotalSeqC",rownames(seurat_object))]
    common_markers=c(common_markers,antibody_name)
  }
  non_common=setdiff(unique(as.character(all.markers$gene)),common_markers)
  all.markers=subset(all.markers,gene%in%non_common)
  top5 <- all.markers %>% group_by(cluster) %>% top_n(5, p_val_adj)
  top20 <- all.markers %>% group_by(cluster) %>% top_n(20, p_val_adj)
  heatgene=unique(as.character(top20$gene))
  
  if(!isTRUE(tryCatch(nrow(.get_assay_data(seurat_object, "scale.data"))>0, error=function(e) FALSE))){
    DefaultAssay(seurat_object)="RNA"
    seurat_object=ScaleData(seurat_object)
  }
  
  ave_cluster=AverageExpression(seurat_object,assays = "RNA",slot = "data",return.seurat = T)
  pdf(paste0("heatmap_",filename,"_top20.pdf"),width =15, height = 40)
  f1=DoHeatmap(ave_cluster, features = heatgene,
               label = TRUE,size=3,assay = "RNA",lines.width = 1,
               group.colors = color_crc,draw.lines = F)+
    scale_fill_gradientn(colors = fill_colors)
  print(f1)
  dev.off()
  
  pdf(paste0("heatmap_",filename,"_width_top20.pdf"),width =16, height = 12)
  print(f1)
  dev.off()
  
  pdf(paste0("heatmap_",filename,"_top5.pdf"),width =15, height = 12)
  f2=DoHeatmap(object = ave_cluster, features = unique(as.character(top5$gene)),
               label = TRUE,size=3,group.colors = color_crc,draw.lines = F)+
    scale_fill_gradientn(colors = fill_colors)
  print(f2)
  dev.off()
  
  pdf(paste0("heatmap_",filename,"_width_top5.pdf"),width =16, height = 12)
  print(f2)
  dev.off()
  
  png(paste0("heatmap_",filename,"_top20.png"),width = 1500, height = 4000)
  print(f1)
  dev.off()
  
  png(paste0("heatmap_",filename,"_top5.png"),width = 750/1.5, height = 600/1.5)
  print(f2)
  dev.off()
  
  pdf(paste0("Dotplot_",filename,"_nolegend.pdf"),width =20, height = 5.18)
  f3=DotPlot(object = seurat_object, features = rev(x = unique(top5$gene)), cols = c("blue","firebrick2"), 
             dot.scale = 5) + RotatedAxis()+NoLegend()
  print(f3)
  dev.off()
  
  pdf(paste0("total_heatmap_",filename,"_top5.pdf"),width =15, height = 12)
  f4=DoHeatmap(object = seurat_object, features = unique(top5$gene),
               label = TRUE,size=3,group.colors = color_crc,draw.lines = T)+
    scale_fill_gradientn(colors = fill_colors)
  print(f4)
  dev.off()
  
  pdf(paste0("total_heatmap_",filename,"_top20.pdf"),width =15, height = 40)
  f5=DoHeatmap(object = seurat_object, features = unique(top20$gene),
               label = TRUE,size=3,group.colors = color_crc,draw.lines = T)+
    scale_fill_gradientn(colors = fill_colors)
  print(f5)
  dev.off()
  
  export=list(top5=unique(as.character(top5$gene)),top20=unique(as.character(top20$gene)))
  return(export)
}


## get seurat group GO dot plot ##
#' get seurat group GO dot plot
#'
#' @export
get_seuratgroupGO = function(seurat,filename){
  message('in current for mouse only,and the seurat ident shall be group info')
  all.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
  top50 <- all.markers %>% group_by(cluster) %>% top_n(50, -p_val_adj)
  
  ids=bitr(top50$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
  top50=merge(top50,ids,by.x='gene',by.y='SYMBOL')
  gcSample=split(top50$ENTREZID, top50$cluster) 
  xx <- compareCluster(gcSample,
                       fun = "enrichGO",
                       OrgDb='org.Mm.eg.db', pvalueCutoff = 0.05
  )
  p <- dotplot(xx)
  p <- p + theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))+labs(title = filename)
  pdf(paste0(filename,'_dotplot.pdf'),width = 6,height = 8)
  print(p)
  dev.off()
  
}


## get the prediction of celltype from liver cell atlas biomarkers ##
#' get the prediction of celltype from liver cell atlas biomarkers
#'
#' @export
get_predcelltype = function(seurat_obj = NULL, excel_fullpath = NULL,species = 'Human'){
  message('species should be Human or Mouse')
  sheet_names <- excel_sheets(excel_fullpath)
  cell_type_markers <- list()
  
  if (species == 'Human') {
    sheet_names = grep(sheet_names,pattern='DEPs|Mouse',value = T,invert = T)
  }else{
    sheet_names = grep(sheet_names,pattern='DEPs|Human',value = T,invert = T)
  }
  for (sheet in sheet_names) {
    df <- read_excel(excel_fullpath, sheet = sheet)
    

    genes <- df[[1]]
    

    genes <- na.omit(genes)

    cell_type_markers[[sheet]] <- genes

    message(paste("already read in ", sheet, " - found ", length(genes), " genes"))
  }
  
  
  valid_markers <- list()
  for (cell_type in names(cell_type_markers)) {
    markers <- cell_type_markers[[cell_type]]

    valid_genes <- markers[markers %in% rownames(seurat_obj)]
    
    if (length(valid_genes) > 0) {
      valid_markers[[cell_type]] <- valid_genes
    }
  }

  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = valid_markers,
    name = names(valid_markers)
  )

  score_columns <- paste0(names(valid_markers), "1")
  cell_pred = grep(colnames(seurat_obj@meta.data), pattern = species)
  max_scores <- apply(seurat_obj@meta.data[, cell_pred], 1, which.max)
  predicted_types <- names(valid_markers)[max_scores]
  seurat_obj@meta.data[, cell_pred] = NULL
  seurat_obj$predicted_cell_type <- predicted_types
  seurat_obj$predicted_cell_type = str_replace(seurat_obj$predicted_cell_type,pattern = 'Human |Mouse ',replacement = '')
  seurat_obj$predicted_cell_type = str_replace(seurat_obj$predicted_cell_type,pattern = ' DEGs',replacement = '')
  return(seurat_obj)
}
