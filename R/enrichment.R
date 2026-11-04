# Part of the getBio package.

## GO and kegg analysis ##
#' GO and kegg analysis
#'
#' @export
myGO_KEGG=function(species="m",up_gene,down_gene,width=7.54,height=5.18,show=20,filename){
  require("org.Mm.eg.db")
  require(DOSE)
  require(clusterProfiler)
  require("org.Hs.eg.db")
  
  species=species
  if(species=="h"|species=="human"){
    organism_db="org.Hs.eg.db"
    organism="hsa"
  }else{
    organism_db="org.Mm.eg.db"
    organism="mmu"
  }
  png.width = 100*width
  png.height = 100*height
  pdf.width = width
  pdf.height = height
  up_geneid = bitr(up_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  up_geneid = up_geneid$ENTREZID
  go_pvalue=0.05
  up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db, 
                     keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n=0
  while(length(up_ego)==0&n<10){
    go_pvalue=go_pvalue+0.05
    up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
    n=n+1
  }
  n=0
  if(length(up_ego)==0){
    print("There were no UP GO terms enriched,please check your data!!")
  }else{
    while (nrow(up_ego)<5&n<10){
      go_pvalue=go_pvalue+0.05
      up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db, 
                         keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }
  down_geneid = bitr(down_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  down_geneid = down_geneid$ENTREZID
  go_pvalue=0.05
  down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db, 
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n=0
  while(length(down_ego)==0&n<10){                                                                                                                                                                                                                 go_pvalue=go_pvalue+0.05
  down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n=n+1
  }
  n=0
  if(length(up_ego)==0){
    print("There were no DOWN GO terms enriched,please check your data!!")
  }else{
    while (nrow(down_ego)<5&n<10){
      go_pvalue=go_pvalue+0.05
      down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db, 
                           keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }
  library(clusterProfiler)
  
  up_ego=clusterProfiler::simplify(up_ego,cutoff = 0.7,by = "p.adjust",select_fun = min)
  
  down_ego=clusterProfiler::simplify(down_ego,cutoff=0.7,by="p.adjust",select_fun=min)
  
  
  write.csv(up_ego,paste0(filename,"_GO_UP_result.csv"))
  write.csv(down_ego,paste0(filename,"GO_DOWN_result.csv"))
  
  if (nrow(up_ego)>1) {
    png(paste0(filename,"_GO_UP_top20_dotplot.png"), width = png.width, height = png.height)
    d1=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"),x='Count')
    print(d1)
    dev.off()
    pdf(paste0(filename,"_GO_UP_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d2=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"))
    print(d2)
    dev.off()
    
  }
  
  if (nrow(down_ego)>1) {
    pdf(paste0(filename,"_GO_down_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d3=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d3)
    dev.off()
    png(paste0(filename,"_GO_down_top20_dotplot.png"), width = png.width, height = png.height)
    d4=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d4)
    dev.off()
  }

}


#' myGO_KEGG2
#'
#' @export
myGO_KEGG2=function(species="m",up_gene,down_gene,png.width=1000,png.height=500,pdf.width=7,pdf.height=5,show=20,filename=NULL){
  #test_package(c("org.Hs.eg.db",'org.Mm.eg.db','graphics','utils','DOSE'))
  
  if(species=="h"|species=="human"){
    organism_db="org.Hs.eg.db"
    organism="hsa"
  }else{
    organism_db="org.Mm.eg.db"
    organism="mmu"
  }
  
  #upgene
  up_geneid <- bitr(up_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  up_geneid <- up_geneid$ENTREZID
  go_pvalue <- 0.05
  up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                     keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n <- 0
  while(length(up_ego)==0 & n < 10){
    go_pvalue <- go_pvalue+0.05
    up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
    n=n+1
  }
  
  n <- 0
  if(length(up_ego)==0){
    print("There were no UP GO terms enriched,please check your data!!")
  }else{
    while (nrow(up_ego)<5&n<10){
      go_pvalue <- go_pvalue+0.05
      up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db,
                         keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }
  
  #downgene
  down_geneid <- bitr(down_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  down_geneid <- down_geneid$ENTREZID
  go_pvalue <- 0.05
  down_ego <- enrichGO(gene = down_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  
  n <- 0
  while(length(down_ego)==0&n<10){                                                                                                                                                                                                                 go_pvalue=go_pvalue+0.05
  down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db,
                       keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  n=n+1
  }
  n <- 0
  if(length(down_ego)==0){
    print("There were no DOWN GO terms enriched,please check your data!!")
  }else{
    while (nrow(down_ego)<5&n<10){
      go_pvalue=go_pvalue+0.05
      down_ego <- enrichGO(gene=down_geneid,OrgDb = organism_db,
                           keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
      n=n+1
    }
  }
  
  up_ego <- simplify(up_ego)
  down_ego <- simplify(down_ego)
  utils::write.csv(up_ego,paste0(filename,"_GO_UP_result.csv"))
  utils::write.csv(down_ego,paste0(filename,"GO_DOWN_result.csv"))
  
  if (nrow(up_ego)>1) {
    png(paste0(filename,"_GO_UP_top20_barplot.png"), width = png.width, height = png.height)
    b1=graphics::barplot(up_ego, showCategory=show,title=paste0("GO_","_UP")) #条形图
    print(b1)
    dev.off()
    png(paste0(filename,"_GO_UP_top20_dotplot.png"), width = png.width, height = png.height)
    d1=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"))
    print(d1)
    dev.off()
    pdf(paste0(filename,"_GO_UP_top20_barplot.pdf"), width = pdf.width, height = pdf.height)
    b2=graphics::barplot(up_ego, showCategory=show,title=paste0("GO_","_UP")) #条形图
    print(b2)
    dev.off()
    pdf(paste0(filename,"_GO_UP_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d2=dotplot(up_ego,showCategory=show,title=paste0("GO_","_UP"))
    print(d2)
    dev.off()
  }
  
  if (nrow(down_ego)>1) {
    pdf(paste0(filename,"_GO_down_top20_barplot.pdf"), width = pdf.width, height = pdf.height)
    b3=graphics::barplot(down_ego, showCategory=show,title=paste0("GO_","_DOWN")) #条形图
    print(b3)
    dev.off()
    pdf(paste0(filename,"_GO_down_top20_dotplot.pdf"), width = pdf.width, height = pdf.height)
    d3=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d3)
    dev.off()
    png(paste0(filename,"_GO_down_top20_barplot.png"), width = png.width, height = png.height)
    b4=graphics::barplot(down_ego, showCategory=show,title=paste0("GO_","_DOWN")) #条形图
    print(b4)
    dev.off()
    png(paste0(filename,"_GO_down_top20_dotplot.png"), width = png.width, height = png.height)
    d4=dotplot(down_ego,showCategory=show,title=paste0("GO_","_DOWN"))
    print(d4)
    dev.off()
  }
  
  
  
  kegg_pvalue <- 0.05
  up_ekk <- enrichKEGG(gene=up_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                       pvalueCutoff=kegg_pvalue)
  n <- 0
  while (nrow(up_ekk) < 5 & n < 10){
    kegg_pvalue <- kegg_pvalue+0.05
    up_ekk <- enrichKEGG(gene=up_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                         pvalueCutoff=kegg_pvalue)
    n = n+1
  }
  up_symbol_ekk <- DOSE::setReadable(up_ekk, OrgDb = organism_db,keyType = "ENTREZID")
  
  kegg_pvalue <- 0.05
  down_ekk <- enrichKEGG(gene=down_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                         pvalueCutoff = kegg_pvalue)
  n <- 0
  while (nrow(down_ekk)<5 & n<10){
    kegg_pvalue <- kegg_pvalue+0.05
    down_ekk <- enrichKEGG(gene=down_geneid,keyType = 'ncbi-geneid',organism = organism, pAdjustMethod = "BH",
                           pvalueCutoff=kegg_pvalue)
    n=n+1
  }
  down_symbol_ekk <- DOSE::setReadable(down_ekk, OrgDb = organism_db,keyType = "ENTREZID")
  utils::write.csv(up_symbol_ekk,paste0(filename,"_KEGG_UP_result.csv"))
  utils::write.csv(down_symbol_ekk,paste0(filename,"_KEGG_DOWN_result.csv"))
  
  if (nrow(up_symbol_ekk)>1) {
    png(paste0(filename,"_KEGG_up_cnetplot.png"), width = png.width, height = png.height)
    k1=cnetplot(up_symbol_ekk, categorySize="pvalue")
    print(k1)
    dev.off()
    pdf(paste0(filename,"_KEGG_up_cnetplot.pdf"), width = pdf.width, height = pdf.height)
    k2=cnetplot(up_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k2)
    dev.off()
  }
  
  if (nrow(down_symbol_ekk)>1) {
    png(paste0(filename,"_KEGG_down_cnetplot.png"), width = png.width, height = png.height)
    k3=cnetplot(down_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k3)
    dev.off()
    pdf(paste0(filename,"_KEGG_down_cnetplot.pdf"), width = pdf.width, height = pdf.height)
    k4=cnetplot(down_symbol_ekk, categorySize="pvalue",showCategory=10)
    print(k4)
    dev.off()
  }
  
  if (nrow(up_ekk)>1) {
    png(paste0(filename,"_KEGG_UP_dotplot.png"), width = png.width, height = png.height)
    k5=dotplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k5)
    dev.off()
    pdf(paste0(filename,"_KEGG_UP_dotplot.pdf"), width = pdf.width, height = pdf.height)
    k6=dotplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k6)
    dev.off()
    png(paste0(filename,"_KEGG_UP_barplot.png"), width = png.width, height = png.height)
    k7=graphics::barplot(up_ekk,showCategory = show,title = paste0("KEGG_","_UP" ))
    print(k7)
    dev.off()
    pdf(paste0(filename,"_KEGG_UP_barplot.pdf"), width = pdf.width, height = pdf.height)
    k8=graphics::barplot(up_ekk,showCategory =show,title = paste0("KEGG_","_UP" ))
    print(k8)
    dev.off()
  }
  
  if (nrow(down_ekk)>1) {
    png(paste0(filename,"_KEGG_DOWN_dotplot.png"), width = png.width, height = png.height)
    k9=dotplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k9)
    dev.off()
    pdf(paste0(filename,"_KEGG_DOWN_dotplot.pdf"), width = pdf.width, height = pdf.height)
    k10=dotplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k10)
    dev.off()
    png(paste0(filename,"_KEGG_DOWN_barplot.png"), width = png.width, height = png.height)
    k11=graphics::barplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k11)
    dev.off()
    pdf(paste0(filename,"_KEGG_DOWN_barplot.pdf"), width = pdf.width, height = pdf.height)
    k12=graphics::barplot(down_ekk,showCategory = show,title = paste0("KEGG_","_DOWN" ))
    print(k12)
    dev.off()
  }
  
  
  
  
}


## get the go and kegg result from nrDEG ##
#' get the go and kegg result from nrDEG
#'
#' @export
get_GoKegg = function(nrDEG, filename, width, height,species = 'm'){
  upgene = nrDEG[nrDEG$g=='UP',c('symbol','t')]
  upgene = upgene[order(upgene$t,decreasing = T),1]%>%as.character()
  downgene = nrDEG[nrDEG$g=='DOWN',c('symbol','t')]
  downgene= downgene[order(downgene$t,decreasing = T),1]%>%as.character()
  if (species == 'human') {
    upgene = str_to_upper(upgene)
    downgene = str_to_upper(downgene)
  }
  filename = filename
  myGO_KEGG(species = species,up_gene = upgene,down_gene = downgene,filename = filename,width = width,height = height)
}


## get certain GO pathway from go.csv file ##
#' get certain GO pathway from go.csv file
#'
#' @export
get_goseletCSV = function(filepath, GO = NULL, n = NULL){
  message("if GO is null, then top 'n' will be choose")
  GO_path = read.csv(filepath,header = T,row.names = 1)
  if (is.null(GO)) {
    GO = head(GO_path[order(GO_path$Count,decreasing = T),],n)%>%row.names()
  }
  Choose_go = GO_path[GO,]
  Choose_go = Choose_go[order(Choose_go$p.adjust,decreasing = F),]
  Choose_go$Description = factor(Choose_go$Description,levels=Choose_go$Description)
  return(Choose_go)
}


## get self defined GO barplot from go.csv
#' get self defined GO barplot from go.csv
#'
#' @export
get_gobarplot = function(GO_up,GO_down,filename = NULL){
  GO_down$group = 'down'
  GO_up$group = 'up'
  merge_GO = rbind( GO_up,GO_down)
  merge_GO = merge_GO[order(merge_GO$group,merge_GO$p.adjust,decreasing = F),]
  merge_GO$Description = factor(merge_GO$Description,levels = rev(merge_GO$Description))
  merge_GO$nudge_y = 0.1
  merge_GO$hjust = 0
  
  write.csv(merge_GO,'merge_GO.csv')
  
  p1 = ggplot(merge_GO, aes(x=Description, y=-log10(p.adjust),fill = group))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = c('lightblue','pink'))+
    coord_flip()+
    theme_classic()+  
    theme(
      text = element_text(size=15),   
      axis.title.y=element_blank(), 
      axis.title.x=element_text(size=15),
      axis.text.y = element_blank(),
      legend.title=element_blank() 
    )+scale_y_continuous(expand = c(0,0))+
    geom_text(aes(label = Description, y = nudge_y),
              nudge_x =0,nudge_y =0,hjust = merge_GO$hjust,
              size = 5)
  print(p1)
  ggsave(p1,width = 6,height = 8,filename = paste0(filename, '_GO_barplot.pdf'))
}


## get go pathway rather than do mygokegg ##
#' get go pathway rather than do mygokegg
#'
#' @export
get_gopathway = function(up_gene,species ='mouse'){
  if(species=="h"|species=="human"){
    organism_db="org.Hs.eg.db"
    organism="hsa"
  }else{
    organism_db="org.Mm.eg.db"
    organism="mmu"
  }
  up_geneid = bitr(up_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism_db)
  up_geneid = up_geneid$ENTREZID
  go_pvalue=0.05
  up_ego <- enrichGO(gene=up_geneid,OrgDb = organism_db, 
                     keyType = 'ENTREZID', ont="BP",pvalueCutoff=go_pvalue,readable=TRUE)
  return(up_ego@result)
}
