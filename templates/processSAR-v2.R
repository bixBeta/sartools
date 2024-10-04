#!/usr/bin/env Rscript

arg <-  commandArgs(trailingOnly = T)

suppressPackageStartupMessages(library(dplyr))


if (length(arg)==0) {
  print(" Usage = Rscript processSAR.R < PIN > ")  
  stop("Please provide the PIN !!! \n", call.=FALSE)
  
} 


library(progress)

pb <- progress_bar$new(total = 100)


for (i in 1:20) {
  pb$tick()
  Sys.sleep(1 / 10)
}




dir = paste0(getwd(), "/")

fileNames <- list.files( dir, pattern = ".complete")
filePath <- paste0(dir, fileNames)





# import SAR tool files as objects 

contrasts.list = list()
for (i in 1:length(fileNames)) {
  
  contrasts.list[[i]] <- read.table(file= filePath[i], header = T, sep = "\t")
  names(contrasts.list)[[i]] <- strsplit(fileNames[i], "\\.")[[1]][1]

}


for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}

z <- c("dispGeneEst","dispFit","dispMAP","dispersion","betaConv","maxCooks" )
contrasts.list.rm.z = lapply(X = contrasts.list, FUN = function(x){
  x %>% select(-all_of(z))
})

for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}


contrasts.list.final = list()

for (i in 1:length(contrasts.list.rm.z)) {
  contrasts.list.final[[i]] <- 
    contrasts.list.rm.z[[i]] %>% rename(!!paste0(names(contrasts.list.rm.z)[[i]],".FoldChange") := FoldChange,
                                   !!paste0(names(contrasts.list.rm.z)[[i]],".log2FoldChange") := log2FoldChange,
                                   !!paste0(names(contrasts.list.rm.z)[[i]],".stat") := stat,
                                   !!paste0(names(contrasts.list.rm.z)[[i]],".pvalue") := pvalue,
                                   !!paste0(names(contrasts.list.rm.z)[[i]],".padj") := padj)
  
  names(contrasts.list.final)[[i]] <- names(contrasts.list.rm.z)[[i]]
}
  
ref.df = contrasts.list.final[[1]]

if (length(contrasts.list.final) == 1) {
  
  
  write.table(ref.df, paste0(arg[1],".final.txt"), sep = "\t", quote = F, row.names = F)
  system(paste("mv *.final.txt ../ "))
  
  for (i in 1:20) {
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  message("")
  message("Only 1 DESeq2 contrast available")
  
  
} else {
   
  tables.for.join = lapply(contrasts.list.final[-1], function(x){
    
    x = x %>% data.frame %>% select(Id,matches("vs"))
    
  })
  
  table.joined = suppressMessages(Reduce(f = full_join, tables.for.join))
  
  final.df = left_join(ref.df, table.joined, by = "Id")
  write.table(final.df, paste0(arg[1],".final.txt"), sep = "\t", quote = F, row.names = F)
  
  system(paste("mv *.final.txt ../ "))
  for (i in 1:20) {
    pb$tick()
    Sys.sleep(1 / 100)
  }
  message("")
  message("More than 1 DESeq2 contrast's available")
  
}


