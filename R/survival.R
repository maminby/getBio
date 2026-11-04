# Part of the getBio package.

##### get the associated genes from choose_gene #####
#' get the associated genes from choose_gene
#'
#' @export
Get_associated_gene = function(exprSet,choose_gene){
  library(dplyr)
  library(tidyr)
  y <- as.numeric(exprSet[as.character(choose_gene),])
  rownames <- rownames(exprSet)
  cor_data_df <- data.frame(rownames)
  for (i in 1:length(rownames)){
    test <- cor.test(as.numeric(exprSet[i,]),y,type="spearman")
    cor_data_df[i,2] <- test$estimate
    cor_data_df[i,3] <- test$p.value
  }
  names(cor_data_df) <- c("symbol","correlation","pvalue")
  
  cor_data_sig <- cor_data_df %>% 
    filter(pvalue < 0.05) %>% 
    arrange(desc(abs(correlation)))%>% 
    dplyr::slice(1:500)
  write.table(cor_data_sig,"Survival_associated_gene.txt",sep = '\t')
  return(cor_data_sig)
}
