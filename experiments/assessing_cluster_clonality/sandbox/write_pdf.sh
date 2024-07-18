#!/bin/bash

pattern=$1

for sample in ${pattern}*.gv;do
    length=${#sample}
    echo "${sample:0:length-3}.pdf"
    dot -Tpdf -o ${sample:0:length-3}.pdf ${sample}
done    
