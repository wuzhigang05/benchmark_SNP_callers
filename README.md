## Benchmark Three SNP callers (Samtools/bcftools, GATK, VariantTools package in R) ##

SNP is the acronym for Single Nucleotide Polymorphism. SNP is associated with human's suspecility to
certain diseases, response to drug and environment pathogens. Therefore finding accurate SNP information 
is very important. However, the fact that many tools exist and that the agreement between 
different SNP calllers are little. In this project, we benchmared three popular SNP 
callers in terms of sensitivity and specifity. 
These three tools are: [samtools/bcftools](http://samtools.sourceforge.net/samtools.shtml), 
[GATK](http://bit.ly/1p10oNM)
and [VariantTools](http://www.bioconductor.org/packages/release/bioc/html/VariantTools.html).

### Requirements ###

This project is meant to run ONLY in a Linux environment. Specifically, you should be able to run it without any 
problem in [UCR biocluster platform](http://manuals.bioinformatics.ucr.edu/home/hpc). If you want to run this under other Linux environment, you have to install below requirements. 

* Java Jars
 * GenomeAnalysisTK.jar version 3.1.1 from [GATK](http://bit.ly/1p10oNM)
 * [PICARD](http://picard.sourceforge.net/) project 
    * AddOrReplaceReadGroups.jar  
    * CreateSequenceDictionary.jar
    * MarkDuplicates.jar
* [samtools version 0.1.19](http://sourceforge.net/projects/samtools/)
* [seqtk](https://github.com/lh3/seqtk)
* python packages:
  * [numpy](http://www.numpy.org/)
  * [pandas](http://pandas.pydata.org/)
  * [matplotlib](http://matplotlib.org/)
* R packages:
 * [VariantTools](http://www.bioconductor.org/packages/release/bioc/html/VariantTools.html)
 * [gmapR](http://www.bioconductor.org/packages/2.12/bioc/html/gmapR.html)
 * [rtracklayer](http://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)

### Usage ###
After git clone this project, just type `make` in your terminal, it will perform all the necessary 
steps (downloading test data, performing alignments, call SNP using different SNP callers, SNP 
fileration and etc) and finally generate a ROC curve pdf named "plot.pdf" in your current direcory.
WARNING, make sure you have 10GB diskspace available before you type ``make`` because of the intermediate fastq files.
```bash
git clone http://github.com/wuzhigang05/benchmark_SNP_callers 
make
make open
```
### Results ###
* overlap among three SNP callers
![alt tag]()

* TPR and FDR of three SNP callers
![alt tag]()

### Comments ###
Send your comments to zhigang dot wu at email dot ucr dot edu

