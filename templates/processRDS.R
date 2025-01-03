#!/usr/bin/env Rscript

args <-  commandArgs(trailingOnly = T)
# check for required argument

if (length(args)==0) {
  print(" Usage = Rscript processRDS.R < SARTools RDATA >")
  stop("Missing SARTools RDATA !!! \n", call.=FALSE)
  
}

suppressPackageStartupMessages(library(SARTools))
suppressPackageStartupMessages(library(dplyr))

pin <- args[1]
load(args[2])

dds <- out.DESeq2$dds
vst = varianceStabilizingTransformation(object = dds, blind = T)
vst

getPCAs= function(vst_, target_){
  
  
  meta = as.data.frame(target_)
  rv <- rowVars(assay(vst_))
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  
  pca <- prcomp(t(assay(vst_)[select,]))
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  pVar.df <- as.data.frame(percentVar)
  pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
  
  pVar.df = pVar.df[ , order(names(pVar.df))]
  pVar.df$percentVar = pVar.df$percentVar * 100
  pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
  
  
  d <- data.frame(pca$x, label=rownames(pca$x))
  d2 <- left_join(d, meta, by = "label")
  
  
  return(list(
    prcomp.out = pca,
    Variance.df    = pVar.df,
    colData    = meta,
    PCA.df      = d2
  ))
  
}

pca_rds = getPCAs(vst_ = vst, target_ = target)

suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
pc1 = ggplot(pca_rds$PCA.df, aes(x=PC1, y=PC2, color = group)) +
  geom_point(size=2.5) +
  geom_label_repel(aes(label = label),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey55', show.legend = F) + 
  xlab(paste0(pca_rds$Variance.df$x[1], "  ", pca_rds$Variance.df$percentVar[1], "%") ) +
  ylab(paste0(pca_rds$Variance.df$x[2], "  ", pca_rds$Variance.df$percentVar[2], "%") )

png(paste0(pin, "_PC1_PC2.png"), width = 1200, height = 1200, res = 150)
pc1
dev.off()


write.csv(pca_rds$PCA.df, file = paste0(pin, "_EIGENVALUES.csv"))
write.table(target, file = paste0(pin, "_targetFile.txt"), quote = F)
saveRDS(pca_rds, file = paste0(pin, "_PCA_EIGEN.RDS"))