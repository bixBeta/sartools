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
params.gbcov            = "${launchDir}/GBCOV"

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
        label "sar_tools"

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

            template 'bash-template.sh' 
            




}



process QMD {

       tag "${id}"
       label "qmds"
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
        label "gbcov"

        input:
            path(gbcov)

        
        output:

            path "*.png"            , emit: gbpng

        when:
            params.quarto == "k1"

        
        script: 

            template "pdf2png.sh"



}





workflow  NOGBC {

        SARTOOLS(params.id, params.ref, ch_target, ch_counts)
        
        ch_figures = SARTOOLS.out.figures
        
        ch_all     = SARTOOLS.out.sartoolsOut
        
        QMD(params.id, params.ref, ch_target, ch_figures, params.quarto, params.genome, params.annots)

}


workflow  GBC {

        SARTOOLS(params.id, params.ref, ch_target, ch_counts)
        GBCOV(params.gbcov)

        ch_figures = SARTOOLS.out.figures
                        .concat(GBCOV.out.gbpng)
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