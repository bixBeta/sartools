#!/bin/bash

echo '""""""""""""""""""""""""""""""""""""""'
echo  ""
echo  "Hello from the bash template script"
echo  "project id is ${id}"
echo '""""""""""""""""""""""""""""""""""""""'



Rscript ${projectDir}/SAR_tools.R ${id} ${ref} ${target} ${launchDir}