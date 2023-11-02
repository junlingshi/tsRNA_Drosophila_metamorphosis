library(AnnotationDbi)
library( org.Dm.eg.db)
library(clusterProfiler)                                

count<-up$Geneid
gene <-  mapIds(org.Dm.eg.db, count, 'ENTREZID', 'FLYBASE')
BP.params <- enrichGO(   gene   = gene,    
OrgDb  = org.Dm.eg.db,    
ont   = "BP"  ,    
pAdjustMethod = "BH",    
pvalueCutoff  = 0.01,    
qvalueCutoff  = 0.05) 
BP.list <- setReadable(BP.params, org.Dm.eg.db, keyType = "ENTREZID")
dotplot(BP.list, showCategory=30)
