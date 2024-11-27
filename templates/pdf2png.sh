#!/bin/bash

cd ${params.gbcov}
convert -density 300 \\*.geneBodyCoverage.curves.pdf curves.png
convert -density 300 \\*.geneBodyCoverage.heatMap.pdf heatMap.png


