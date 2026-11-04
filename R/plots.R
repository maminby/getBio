# Part of the getBio package.

## get correlation among samples ##
#' get correlation among samples
#'
#' @export
get_correlation = function(data){
  require('RColorBrewer')
  require('pheatmap')
  require("ggplot2")
  require('limma')
  require('psych')
  require('ggrepel')
  require('scatterplot3d')
  pearson=cor(data)
  pearson_t=corr.test(pearson, adjust = "none", use = "complete")
  pearson_matrix=pearson_t$r
  colors <- colorRampPalette( rev(brewer.pal(9, "RdYlBu")))(255)
  pdf("samples_correlation.pdf", width =7.54, height = 5.18)
  pheatmap(pearson_matrix,
           cluster_cols = TRUE,cluster_rows = TRUE,
           color = colors,treeheight_row = 10,treeheight_col = 10,
           clustering_distance_row = "correlation",
           clustering_distance_cols = "correlation",
           main= "correlation",fontsize_col = 5,fontsize_row = 5)
  dev.off()
}


## get PCA among samples ##
#' get PCA among samples
#'
#' @export
get_pca = function(data,group_list){
  require('RColorBrewer')
  require('pheatmap')
  require("ggplot2")
  require('limma')
  require('psych')
  require('ggrepel')
  require('scatterplot3d')
  group_list=as.character(unlist(group_list))
  data=data[rowSums(data<=0.1)<ncol(data),]
  data = t(data)
  pc <- prcomp(data,scale=T)
  pcx <- data.frame(pc$x)
  pcr <- cbind(samples = rownames(pcx),group_list,pcx) 
  p <- ggplot(pcr, aes(PC1, PC2)) + geom_point(size = 3,aes(color=group_list))+theme_bw()+theme()+
    theme(panel.grid=element_blank(),axis.line=element_line(size=0.1,colour="black"))+
    scale_color_manual(values = color_crc)
  #geom_text_repel(aes(label = samples),hjust = 0.5, vjust = 2,size=3)
  pdf("pca_2D.pdf", width =7.54, height = 5.18)
  print(p)
  dev.off()
  p1 <- p + stat_ellipse(level = 0.95, show.legend = F)
  pdf("pca_2D_ellipse.pdf", width =7.54, height = 5.18)
  print(p1)
  dev.off()
}


## volcano plot ##


#' get_volcano
#'
#' @export
get_volcano <- function(nrDEG, genes_list = NULL, save = F, filename = NULL) {
  library(ggplot2)
  library(ggrepel)
  p <- ggplot(nrDEG, aes(x = logFC, y = AveExpr)) +
    geom_point(data = subset(nrDEG, adj.P.Val >= 0.05), aes(size = 0.1, color = "Non-significant"), alpha = 0.6) +
    geom_point(data = subset(nrDEG, adj.P.Val < 0.05), aes(size = -log10(adj.P.Val), color = ifelse(logFC > 1, "Upregulated", ifelse(logFC < -1, "Downregulated", "Non-significant"))), alpha = 0.6) +
    scale_size_continuous(range = c(0.1, 10)) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "grey")) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.1),panel.grid = element_blank(),panel.background = element_blank())+
    labs(title = paste0(filename,'Volcano Plot'),
         x = "logFC",
         y = "AveExpr",
         size = "-log10(adj.P.Val)",
         color = "Significance")
  
  genes_list <- intersect(nrDEG$symbol,genes_list)
  if (!is.null(genes_list)) {
    nrDEG$annotate <- ifelse(nrDEG$symbol%in%genes_list, nrDEG$symbol, NA)
    p <- p + geom_text_repel(data = nrDEG,aes(label = annotate), 
                             box.padding = unit(0.3, "lines"),
                             point.padding = unit(0.4, "lines"),
                             show.legend = F, size = 4)
  }
  up_regulated <- sum(nrDEG$logFC > 1 & nrDEG$adj.P.Val < 0.05)
  down_regulated <- sum(nrDEG$logFC < -1 & nrDEG$adj.P.Val < 0.05)
  p+labs(color = paste("Significance:","\nUpregulated:", up_regulated, "\nDownregulated:", down_regulated))
  return(p)
  if (save) {
    ggsave(paste0(filename,'_volcano.pdf'))
  }
}


## get certain gene's expression profile in violin plot##
#' get certain gene's expression profile in violin plot
#'
#' @export
get_violinplot = function(data, gene, group_list = group_list,pairinfo = c(),method = "wilcox.test",save = F){
  library(ggpubr)
  library(ggplot2)
  data = cbind(data.frame(t(data),stringsAsFactors = F),group = group_list)
  if (method == "wilcox.test") {
    pairinfo = pairinfo
    for (i in 1:179) {
      pairinfo = c(pairinfo,rep(i,2))
    }
    pairinfo = as.factor(pairinfo)
    p <- ggplot(data, aes_string(x = 'group',y = gene,fill='group_list'))+
      geom_violin()+
      geom_point(size=2, alpha=0.5)+
      geom_line(aes(group=pairinfo), colour="black", linetype="11")+
      xlab("")+
      ylab(paste("Expression of ",gene))+
      theme_classic()+
      scale_y_continuous(expand = c(0,0))+
      theme(legend.position = "none")+
      stat_compare_means(method = "wilcox.test")
    print(p)
  }else{
    p <- ggplot(data, aes_string(x = 'group',y = gene,fill='group_list'))+
      geom_violin()+
      geom_point(size=2, alpha=0.5)+
      scale_y_continuous(expand = c(0,0))+
      geom_line(aes(group=pairinfo), colour="black", linetype="11")+
      xlab("")+
      ylab(paste("Expression of ",gene))+
      theme_classic()+
      theme(legend.position = "none")+
      stat_compare_means(method = "t.test")
    print(p)
  }
  if (save == T) {
    pdf(paste0(gene,'.pdf'), width = 8,height = 8)
    print(p)
    dev.off()
  }
}

#get_survival_bygene = function(expdata,pdata,gene,save = F){
#  library(survival)
#  library(stringr)
#  library(survminer)
#  group_list = ifelse(Get_logic(expdata[gene,]),paste0(gene,'hi'),paste0(gene,'low'))%>%as.character()
#  group_list = data.frame(group_list,row.names = colnames(expdata))
#  pdata <<- cbind(pdata,group_list)
#  pdata$time <<- as.numeric(pdata$survival.time.months.)
 # pdata$Status <<- as.numeric(as.factor(pdata$death.at.fu))
 # fit<- survfit(Surv(time, Status) ~ group_list, data = pdata)
  #
#  p = ggsurvplot(fit, pval = TRUE,  
#                 conf.int = F)
#  print(p)
#  if (save == T) {
#    pdf(paste0(gene,'survival.pdf'),width = 10,height = 10)
#    print(p)
#    dev.off()
#  }
#}


## get the correction of two genes ##
#' get the correction of two genes
#'
#' @export
Get_Scatterplot2 = function(mat.object=NULL,seurat_object=NULL,log=F, genes=c('x','y'),filename=NULL,type=NULL){
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
          geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0,method = type) + geom_rug(color='#7fc97e') +
          theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
          theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
          theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
          xlab(paste(xlab)) + ylab(paste(ylab))
      }else{
        p <- ggplot(data = list1,aes(x = x , y = y)) +
          geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0) + geom_rug(color='#7fc97e') +
          theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
          theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
          theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
          xlab(paste(xlab)) + ylab(paste(ylab))
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
      geom_point(color='#9e4aa8') + geom_smooth(color='#fbc08b',size=2.0,method = 'lm') + geom_rug(color='#7fc97e') +
      theme_bw() +theme(axis.text.x = element_text(face = 'bold',size = 12), axis.text.y = element_text(face = 'bold',size = 12)) +
      theme(axis.title.x = element_text(face = 'bold',size = 12),axis.title.y = element_text(face = 'bold',size = 12)) +
      theme(panel.border=element_blank()) +theme(axis.ticks=element_blank())+
      xlab(paste(xlab)) + ylab(paste(ylab))
    
    pdf(paste0("scatterplot_",filename,".pdf"),width =4, height = 4)
    print(p)
    dev.off()
    
    p <- p + stat_cor(data=list1, method = "spearman")
    pdf(paste0("scatterplot_spearman",filename,".pdf"),width =4, height = 4)
    print(p)
    dev.off()
  }
  
}


## get venn plot ##
#' get venn plot
#'
#' @export
get_vennplot = function(genelist,filename){
  library(eulerr)
  library(VennDiagram)
  if (is.null(names(genelist))) names(genelist) = LETTERS[seq_along(genelist)]
  vennlist = get.venn.partitions(genelist)
  # build the region labels (e.g. "A&B") from the logical membership columns,
  # instead of string-cleaning the ..set.. column
  setnames = names(genelist)
  lab = apply(vennlist[, setnames, drop = FALSE], 1,
              function(r) paste(setnames[as.logical(r)], collapse = "&"))
  vennlist = vennlist[lab != "", , drop = FALSE]
  vennlist$..set.. <- lab[lab != ""]
  euler.obj = euler(setNames(vennlist$..count..,vennlist$..set..))
  
  p1 = plot(euler.obj,
           fills = list(fill=c("red","blue",
                                    "green","darkgreen",
                                    "pink","black",
                                    "purple"),
                                    alpha=0.5),quantities = euler.obj$original.values,
           edges = list(col="white",alpha=0))
  ggsave(p1,filename = paste0(filename,'_venn_plot.pdf'))
  p2 = plot(euler.obj,
            fills = list(fill=c("red","blue",
                                     "green","darkgreen",
                                     "pink","black",
                                     "purple"),
                                     alpha=0.5),
            edges = list(col="white",alpha=0),label = list(col="white",alpha=0))
  ggsave(p2,filename = paste0(filename,'_notext_venn_plot.pdf'))
  print(p1)
  return(vennlist)
}
