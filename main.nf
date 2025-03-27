nextflow.enable.dsl=2


// Global Params: controlled by --prefix 

params.target           = "targetFile.txt"
params.id               = "TREX_id"
params.genome           = "null"
params.annots           = "ENSEMBL"       
params.ref              = null 
params.quarto           = "k3"
params.counts           = "${launchDir}/rawCounts/*rawCounts"
params.help             = false
params.gbcov            = "${launchDir}/GBCOV/*pdf"

// Input Channels: 

ch_counts = channel.fromPath(params.counts) 
                | collect
                

ch_target = channel.fromPath(params.target)


meta_ch = ch_target
                |  splitCsv( header:true , sep:"\t" )
                |  map { row -> [row.label, [row.files, row.group]] }
                |  view


if( params.help ) {

log.info """
S A R  T O O L S     W O R K F L O W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/sartools -r main < args ... >

Args:
    * --id             : TREx Project ID 
    * --ref            : Base Level (Denominator for log2FC calcs, must be in the group column of targetFile)
    * --target         : targetFile.txt (tab delim file with label, files and group mandatory columns)
    * --genome         : Reference genome (GRCh38, GRCm38 etc.)
    * --quarto         : < default: k3 > (render params used in makeReport.sh)
    * --annots         : Annotations source (NCBI or < default: ENSEMBL >)

"""

    exit 0
}

process SARTOOLS {

        tag "${id}"
        label "process_sartools"

        publishDir "${id}_SAR-Tools/"       , mode: 'symlink', overwrite: true, pattern: "*RData"
        publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*_rawCounts.txt"

        input:
            val(id)
            val(ref)
            path(target)
            path(ch_counts)

        
        output:
            // path "*"                              , emit: sartoolsOut
            path "figures/*png"                   , emit: figures
            path "tables/*txt"                    , emit: tables
            path "*_rawCounts.txt"                , emit: rawmatrix
            path "*RData"                         , emit: sarimage

        script:

            // template 'bash-template.sh' 

        """

        #!/usr/bin/env Rscript
            #args <-  commandArgs(trailingOnly = T)
            pin <- "${id}"
            ref <- "${ref}"
            targetFile <- "${launchDir}/targetFile.txt"
            wd <- "${launchDir}"
            templates <- "${projectDir}"

            .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

            
            library(SARTools)
            getwd()
            workDir <- setwd(paste0(wd, "/rawCounts"))      								# working directory for the R session
            
            projectName <- pin                      				# name of the project
            author <- "RSC"                                			# author of the statistical analysis/report
            
            # targetFile <- "targetFile.txt"                           # path to the design/target file
            rawDir <- getwd()                                      # path to the directory containing raw counts files
            #featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
            #                     "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
            #                    "not_aligned", "too_low_aQual")# NULL if no feature to remove
            
            featuresToRemove <- NULL
            varInt <- "group"                                    # factor of interest
            condRef <- ref                                    # reference biological condition
            
            
            
            
            batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example
            
            fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
            cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
            independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
            alpha <- 0.05                                        # threshold of statistical significance
            pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"
            
            typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
            locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors
            
            
            
            #colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            #            "MediumVioletRed","SpringGreen", c2)
            
            
            
            # loading target file
            target <- loadTargetFile(targetFile= targetFile, varInt=varInt, condRef=condRef, batch=batch)
            
            library("viridisLite")
            #colors <- viridisLite::viridis(n = nrow(target), option = "inferno", begin = 0, end = 1 )
            
            colors<-  c("#EF8A62",
            	     "#1f78b4",
            	     "#1b9e77",
            	     "purple3",
            	     "khaki4",
            	     "#E9A3C9",
            	     "#A1D76A",
            	     "#FFFF33",
            	     "grey",
            	     "#b3e2cd",
            	     "#67A9CF",
            	     "peachpuff2",
            	     "red",
            	     "magenta3",
            	     "blue",
            	     "lightyellow",
            	     "black",
            	     "lightblue",
            	     "#E8AEB7",
            	     "#B8E1FF",
            	     "#A9FFF7",
            	     "#94FBAB",
            	     "#82ABA1",
            	     "#1A181B",
            	     "#564D65"
            )
            
            #colors <- c(c1,c2)
            forceCairoGraph <- FALSE
            
            ################################################################################
            ###                             running script                               ###
            ################################################################################
            setwd(workDir)
            library(SARTools)
            if (forceCairoGraph) options(bitmapType="cairo")
            
            # checking parameters
            checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                                   rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                                   condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                                   independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                                   typeTrans=typeTrans,locfunc=locfunc,colors=colors)
            
            # loading target file
            # target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
            
            # loading counts
            counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
            
            # description plots
            majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)
            
            # analysis with DESeq2
            out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                                     locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                                     cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
            
            # PCA + clustering
            exploreCounts(object=out.DESeq2\$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)
            
            # summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
            summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                                      independentFiltering=independentFiltering,
                                                      cooksCutoff=cooksCutoff, alpha=alpha)
            
            # save image of the R session
            save.image(file=paste0(projectName, ".RData"))
            
            # generating HTML report
            writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                               majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                               targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                               condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                               independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                               typeTrans=typeTrans, locfunc=locfunc, colors=colors)
            
            
            
            #load()
            write.table(counts, paste0(projectName, "_rawCounts.txt"), sep = "\t", quote = F, col.names = NA)
            ################################################

        """

}

process VS {

        tag "${id}"
        label "process_sartools"

        // publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*.txt"
        // publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*_EIGENVALUES.csv"
        // publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*_PC1_PC2.png"

        input:
            val(id)
            path(tables)

            // val(ref)
            // path(target)
            // path(ch_counts)

        output:

            path "*txt" , emit:txts


        script:

            """

            
            mkdir others
            mv  *.txt others
            cd others
            mv *complete* ..
            cd ..

            for i in *complete*; do mv \$i `echo \$i | sed 's/vs/_vs_/g'` ; done

            """

}

process RDS {

        tag "${id}"
        label "process_sartools"

        publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*.txt"
        publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*_EIGENVALUES.csv"
        publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*_PC1_PC2.png"

        input:
            val(id)
            path(rdata)

        output:

            path "*.csv"     , emit: eigen
            path "*.RDS"     , emit: eigenrds
            path "*.png"     , emit: ggpca    


        script:

            """
                #!/usr/bin/env Rscript
            
                .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

                suppressPackageStartupMessages(library(SARTools))
                suppressPackageStartupMessages(library(dplyr))

                pin <- "${id}"
                load("${rdata}")

                dds <- out.DESeq2\$dds
                vst = varianceStabilizingTransformation(object = dds, blind = T)
                vst

                getPCAs= function(vst_, target_){
                
                
                meta = as.data.frame(target_)
                rv <- rowVars(assay(vst_))
                select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
                
                pca <- prcomp(t(assay(vst_)[select,]))
                
                percentVar <- pca\$sdev^2 / sum( pca\$sdev^2 )
                pVar.df <- as.data.frame(percentVar)
                pVar.df\$x = as.factor(paste0("PC",rownames(pVar.df)))
                
                pVar.df = pVar.df[ , order(names(pVar.df))]
                pVar.df\$percentVar = pVar.df\$percentVar * 100
                pVar.df\$percentVar = round(pVar.df\$percentVar, digits = 2)
                
                
                d <- data.frame(pca\$x, label=rownames(pca\$x))
                d2 <- left_join(d, meta, by = "label")
                
                
                return(list(
                    prcomp.out = pca,
                    Variance.df    = pVar.df,
                    colData    = meta,
                    PCA.df      = d2
                ))
                
                }

                pca_rds = getPCAs(vst_ = vst, target_ = target)
                pca_rds\$rawCounts = counts(dds, normalized = F)

                suppressPackageStartupMessages(library(ggrepel))
                suppressPackageStartupMessages(library(ggplot2))
                pc1 = ggplot(pca_rds\$PCA.df, aes(x=PC1, y=PC2, color = group)) +
                geom_point(size=2.5) +
                geom_label_repel(aes(label = label),
                                box.padding   = 0.35, 
                                point.padding = 0.5,
                                segment.color = 'grey55', show.legend = F) + 
                xlab(paste0(pca_rds\$Variance.df\$x[1], "  ", pca_rds\$Variance.df\$percentVar[1], "%") ) +
                ylab(paste0(pca_rds\$Variance.df\$x[2], "  ", pca_rds\$Variance.df\$percentVar[2], "%") )

                png(paste0(pin, "_PC1_PC2.png"), width = 1200, height = 1200, res = 150)
                pc1
                dev.off()


                write.csv(pca_rds\$PCA.df, file = paste0(pin, "_EIGENVALUES.csv"))
                write.table(target, file = paste0(pin, "_targetFile.txt"), quote = F, sep = "\t", row.names = F)
                saveRDS(pca_rds, file = paste0(pin, "_PCA_EIGEN.RDS"))



            """

}

process SAR {

        tag "${id}"
        label "process_sartools"

        publishDir "${id}_SAR-Tools/Reports", mode: 'symlink', overwrite: true, pattern: "*.txt"


        input:
            val(id)
            path(txts)

            // val(ref)
            // path(target)
            // path(ch_counts)

        output:

            path "*txt" , emit: finalsar



        script:

            """

            #!/usr/bin/env Rscript

            .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")


            suppressPackageStartupMessages(library(dplyr))



            dir = paste0(getwd(), "/")

            fileNames <- list.files( dir, pattern = ".complete")
            filePath <- paste0(dir, fileNames)

            # import SAR tool files as objects 

            contrasts.list = list()
            for (i in 1:length(fileNames)) {
            
            contrasts.list[[i]] <- read.table(file= filePath[i], header = T, sep = "\t")
            names(contrasts.list)[[i]] <- strsplit(fileNames[i], '\\\.')[[1]][1]

            }


            z <- c("dispGeneEst","dispFit","dispMAP","dispersion","betaConv","maxCooks" )
            contrasts.list.rm.z = lapply(X = contrasts.list, FUN = function(x){
            x %>% select(-all_of(z))
            })



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
            

            
            message("")
            message("Only 1 DESeq2 contrast available")
            
            
            } else {
            
            tables.for.join = lapply(contrasts.list.final[-1], function(x){
                
                x = x %>% data.frame %>% select(Id,matches("vs"))
                
            })
            
            table.joined = suppressMessages(Reduce(f = full_join, tables.for.join))
            
            final.df = left_join(ref.df, table.joined, by = "Id")
            write.table(final.df, paste0(arg[1],".final.txt"), sep = "\t", quote = F, row.names = F)
            
            # system(paste("mv *.final.txt ../ "))

            message("")
            message("More than 1 DESeq2 contrast's available")
            
            }




            """

}



process QMD {

       tag "${id}"
       label "process_quarto"
       publishDir "${id}_SAR-Tools/Reports", mode: 'copy', pattern: "*html"

       input:
            val(id)
            val(ref)
            path(target)
            path(figs)
            val(quarto)
            val(genome)
            val(annots)
            val(qmd)

    output:

        path "*html"                , emit: qmdown

    
    script:

        """
        export XDG_CACHE_HOME=/tmp/quarto-cache
        mkdir -p /tmp/quarto-cache
        cat ${qmd} > test.qmd
        quarto render test.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

        """



}

process GBCOV {

        tag "${id}"
        label "process_sartools"

        input:
            path(pdf)

        
        output:

            path "*png"            , emit: gbpng

        when:
            params.quarto == "k1"

        
        script: 

            """
            convert -density 300 ${pdf[0]} curves.png
            convert -density 300 ${pdf[1]} heatMap.png

            """



}





workflow  NOGBC {

        SARTOOLS(params.id, params.ref, ch_target, ch_counts)
        
        ch_figures = SARTOOLS.out.figures
        
//        ch_all     = SARTOOLS.out.sartoolsOut
        ch_qmd     = channel.value("${projectDir}/assets/qmds/nogbcov-report-nf.qmd")
        ch_qmd
            | view
        SARTOOLS.out.tables.view()    
       // QMD(params.id, params.ref, ch_target, ch_figures, params.quarto, params.genome, params.annots, ch_qmd)
        VS(params.id, SARTOOLS.out.tables)
        VS.out.txts.view()

        RDS(params.id, SARTOOLS.out.sarimage)
        SAR(params.id, VS.out.txts)
}


workflow  GBC {

        SARTOOLS(params.id, params.ref, ch_target, ch_counts)

        ch_pdf = channel.fromPath(params.gbcov)
                    | collect
                    | view

        GBCOV(ch_pdf)

        ch_figures = SARTOOLS.out.figures
                        .concat(GBCOV.out.gbpng)
                        .collect()
                        .view()
        
        ch_all     = SARTOOLS.out.sartoolsOut
        
        QMD(params.id, params.ref, ch_target, ch_figures, params.quarto, params.genome, params.annots)

}


workflow {

    if ( params.quarto == "k1") {

        GBC()
    } 

    else {

        NOGBC()
    }
}
