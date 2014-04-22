#!/bin/bash


function filterTag () {
  if [[ $# -ne 2 ]]
  then
    echo "Usage: ${FUNCNAME} <vcfFile> <tagName> <cutoff>"
    echo "e.g. ${FUNCNAME} sam.raw.vcf DP 3"
    echo "while will only retain all rows with DP >=3"
  exit ${ERROR}
  fi
    
}
