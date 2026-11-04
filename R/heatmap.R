# Part of the getBio package.

## get group heatmap ##
#' get group heatmap
#'
#' @export
get_groupheatmap = function(tpms,geneset1,group_list,filename,scale = 'row',width = 16,
                            height = 16,cluster_rows = T, cluster_cols =T,annotation_col = NULL,
                            cellwidth = 20,cellheight =10){
  #tpms is a expression matrix
  #geneset1 is the gene you want to show, it can also be used for GSEA/GSVA calucation
  #group_list is the group list
  library(pheatmap)
  library(grDevices)
  library(export)
  if (!is.null(annotation_col)) {
    annotation_col = factor(merge[,2])
  }
  genes1 = tpms[row.names(tpms)%in%geneset1,]%>%as.data.frame()
  genes1 = genes1[rowSums(genes1<0.1)<ncol(genes1),]
  name = unique(group_list)
  new.genes1 = data.frame(array(dim = c(nrow(genes1),length(name))))
  colnames(new.genes1) = unique(group_list)
  row.names(new.genes1) = row.names(genes1)
  for (i in 1:length(name)) {
    new.genes1[,i] = rowMeans(genes1[,group_list==name[i]])
  }
  p =pheatmap(new.genes1, scale = scale,annotation_colors = anno_name,
              color =colorRampPalette(c("blue", "white", "red"))(20),
              cellwidth = cellwidth,cellheight = cellheight,annotation_col = annotation_col,
              cluster_rows = cluster_rows, cluster_cols = cluster_cols,main = filename,
              treeheight_row = 0,treeheight_col = 0)
  print(p)
  ggsave(filename = filename, height = height,width = width)
}


## get group heatmap ##
#' get group heatmap
#'
#' @export
get_heatmap = function(tpms,geneset1,group_list,colbygroup=F,
                       filename,scale = 'row',width = 16,gaps_row = NULL,
                       gaps_col = NULL,height = 16,breaks = NA,
                       cluster_rows = T, cluster_cols =T,ann_colors = NULL,legend_breaks = NA){
  #tpms is a expression matrix
  #geneset1 is the gene you want to show, it can also be used for GSEA/GSVA calucation
  #group_list is the group list
  library(pheatmap)
  library(export)
  library(grDevices)
  genes1 = tpms[row.names(tpms)%in%geneset1,]%>%as.data.frame()
  genes1 = genes1[rowSums(genes1<0.1)<ncol(genes1),]
  name = unique(group_list)
  anno_col = data.frame('Group' = group_list)
  rownames(anno_col) = colnames(genes1)
  if (colbygroup==T) {
    genes1 = genes1[,order(group_list,decreasing = T)]
  }
  p =pheatmap(genes1, scale = scale, annotation_col = anno_col,
              color =colorRampPalette(c("#1514FD", "#FAF9FF", "#FE0101"))(20),
              cluster_rows = cluster_rows, cluster_cols = cluster_cols,
              treeheight_row = 0,treeheight_col = 0,gaps_row = gaps_row,gaps_col = gaps_col,
              annotation_colors = ann_colors,legend_breaks = legend_breaks,breaks = breaks)
  graph2ppt(x=p, file = filename, height = height,width = width, aspectr=4)
}
