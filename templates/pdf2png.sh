#!/bin/bash

# for i in *.pdf 
# do 
# 	iSUB=`echo $i | cut -d "." -f3` 
# 	convert -density 300 $i ${iSUB}.png

# done



convert -density 300 *.geneBodyCoverage.curves.pdf curves.png
convert -density 300 *.geneBodyCoverage.heatMap.pdf heatMap.png


