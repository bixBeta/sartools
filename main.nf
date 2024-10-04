nextflow.enable.dsl=2


// Global Params: controlled by --prefix 

params.target           = "targetFile.txt"
params.id               = "TREX_id"
params.genome           = "null"
params.annots           = "ENSEMBL"       
params.ref              = null 
params.quarto           = "-k3"
params.counts           = "${launchDir}/rawCounts/*rawCounts"


// Input Channels: 

ch_counts = channel.fromPath(params.counts) 
                | collect
                

ch_target = channel.fromPath(params.target)


meta_ch = ch_target
                |  splitCsv( header:true , sep:"\t" )
                |  map { row -> [row.label, [row.files, row.group]] }
                |  view




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







workflow  {

        SARTOOLS(params.id, params.ref, ch_target, ch_counts)
        
        ch_figures = SARTOOLS.out.figures
        ch_all     = SARTOOLS.out.sartoolsOut
        
        QMD(params.id, params.ref, ch_target, ch_figures, params.quarto, params.genome, params.annots)

}