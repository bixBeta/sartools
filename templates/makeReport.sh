#!/bin/bash

display_usage(){
  echo "------------------------------------------------------------------------------------------------------------------"
  echo "run the script using the following syntax:"
  echo "    bash scriptName <-k3> <Report_Title> <Genome> <Annot>"
  echo ""
  echo " -k1  = knit with all headers (including MA-plot)"
  echo " -k2  = knit with all headers Filtered (2 projects e.g. PCA-flt.png)"
  echo " -k3  = knit w/o GeneBodyCov"
  echo " -k4  = knit w/o GeneBodyCov Filtered (2 projects e.g. PCA-flt.png)"
  echo " -a1  = knit atac"
  echo " -a2  = knit atac Filtered (2 projects e.g. PCA-flt.png)"
  echo " -s1  = knit smRNA"
  echo "------------------------------------------------------------------------------------------------------------------"
}



knit_full_report(){

  scp ${launchDir}/templates/qmds/full-report-nf.qmd .
  
  quarto render full-report-nf.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}

knit_full_filtered(){

  scp ${launchDir}/templates/qmds/full-report-nf-filtered.qmd .
  
  quarto render full-report-nf-filtered.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}


knit_nogbcov(){

  scp ${launchDir}/templates/qmds/nogbcov-report-nf.qmd .
  
  quarto render nogbcov-report-nf.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}


knit_nogbcov_filtered(){

  scp ${launchDir}/templates/qmds/nogbcov-filtered-nf.qmd .
  
  quarto render nogbcov-filtered-nf.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}


knit_atac(){

  scp ${launchDir}/templates/qmds/atac-nf.qmd .
  
  quarto render atac-nf.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}


knit_atac_filtered(){

  scp ${launchDir}/templates/qmds/atac-filtered-nf.qmd .
  
  quarto render atac-filtered-nf.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}


knit_smrna(){

  scp ${launchDir}/templates/qmds/smrna.qmd .
  
  quarto render smrna.qmd -P title:${id} -P genome:${genome} -P annot:${annots} -o ${id}-Report.html

  rm *.qmd

}





case ${quarto} in
    -h|--help)
      display_usage
      ;;
    -k1|--knit1)
      knit_full_report
      ;;
    -k2|--knit2)
      knit_full_filtered
      ;;
    -k3|--knit3)
      knit_nogbcov
      ;;
    -k4|--knit4)
      knit_nogbcov_filtered
      ;;

    -a1)
    knit_atac
     ;;

    -a2)
      knit_atac_filtered
      ;;

    -s1)
      knit_smrna
      ;;

     *)
      display_usage
      ;;
esac