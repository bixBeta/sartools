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

        publishDir "${id}_SAR-Tools/Reports", mode: 'copy', pattern: "*.txt"
        publishDir "${id}_SAR-Tools/Reports", mode: 'copy', pattern: "*_EIGENVALUES.csv"
        publishDir "${id}_SAR-Tools/Reports", mode: 'copy', pattern: "*_PC1_PC2.png"

        input:
            val(id)
            val(ref)
            path(target)
            path(ch_counts)

        
        output:
             path "*"                              , emit: sartoolsOut
             path "figures/*png"                   , emit: figures


        script:

            // template 'bash-template.sh' 

        """

        #!/usr/bin/env Rscript
            #args <-  commandArgs(trailingOnly = T)
            pin <- ${id}
            ref <- ${ref}
            targetFile <- ${target}
            wd <- ${launchDir}
            templates <- ${projectDir}
            
            if (length(args)<=1) {
              print(" Usage = Rscript test.R <pin> <base-line>")
              stop("Both arguments must be supplied!!! \n", call.=FALSE)
            
            }
            
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
            
            
            
            ################################################
            
            system(paste0(templates, "/templates/generateRaw.R ", projectName, " *.RData"))
            
            #system(paste("mkdir",projectName ))
            #system(paste("mv figures *.html *.RData tables *.txt", projectName))
            
            setwd("./tables")
            #system(paste("pwd"))
            system(paste0(templates, "/templates/vs2_vs_.sh"))
            system(paste0(templates, "/templates/processSAR-v2.R", " ", projectName))
            
            #setwd(projectName)
            #system(paste0("/Users/faraz/macpro/bin/generateRaw.R ../ ", projectName, " ../*.RData"))
            
            
            setwd("../")
            
            system(paste0(templates, "/templates/processRDS.R", " ", projectName, " " ," *.RData"))

        """




}



process QMD {

       tag "${id}"
       label "process_sartools"
       publishDir "${id}_SAR-Tools/Reports", mode: 'copy', pattern: "*html"

       input:
            val(id)
            val(ref)
            path(target)
            path(figs)
            val(quarto)
            val(genome)
            val(annots)


    output:

        path "*html"                , emit: qmdown

    
    script:

        template "makeReport.sh"



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
        
        ch_all     = SARTOOLS.out.sartoolsOut
        
        QMD(params.id, params.ref, ch_target, ch_figures, params.quarto, params.genome, params.annots)

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
