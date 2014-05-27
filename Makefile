all: plot.pdf
ref = Chr4.fasta
refIndex = $(addsuffix .sa, $(ref))
readFile = sreads.fq
# preads is the file name contains reads simulated from pseudo genome
preads = preads.fq
# rreads is the file name contains reads simulated from reference genome
rreads = rreads.fq
sortMem = 5G
mutation = 1000
altCov = 10
refCov = 0
javaArgs = -Xmx4g

export PATH :=$(shell pwd):/rhome/zwu/bin/:/opt/samtools/0.1.19/bin:/opt/java/jdk1.7.0_17/bin/:$(PATH)

sam.raw.vcf: aln.sorted.rmdup.bam $(ref) 
	samtools mpileup -Iuf  $(ref) $< | bcftools view -Nvcg - > $@

gatk.raw.vcf: realigned.bam 
	java $(javaArgs) -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R $(ref) -I $< -o $@

varianttools.raw.vcf\
varianttools.raw.1.vcf varianttools.raw.2.vcf varianttools.raw.3.vcf \
varianttools.raw.4.vcf varianttools.raw.5.vcf varianttools.raw.6.vcf \
varianttools.raw.7.vcf varianttools.raw.8.vcf varianttools.raw.9.vcf: \
run.R data/Chr4.fasta data/Chr4.dict data/Chr4.fasta.fai aln.sorted.rmdup.bam.bai
	Rscript run.R

aln.sorted.rmdup.bam: aln.sorted.bam
	samtools rmdup -s $< $@

aln.sorted.rmdup.bam.bai: aln.sorted.rmdup.bam
	samtools index $<

aln.sorted.bam: $(refIndex) $(readFile)
	bwa mem $(ref) $(readFile) | samtools view -bSF 0x4 - | samtools sort -m$(sortMem) - aln.sorted
	samtools index $@

${refIndex}: ${ref}
	bwa index $<

pChr4.fasta: pseudogenome.py
	./$< ${ref} 

$(preads): pChr4.fasta
	art_illumina -rs 100 -l 55 -f $(altCov) -i $< -o $(basename $(preads))

$(rreads): Chr4.fasta
	touch $(rreads)

$(readFile): $(preads) $(rreads)
	cat $(preads) $(rreads) > $@
	rm -f $(basename $(rreads))* $(basename $(preads))*

validated.aln.sorted.bam: aln.sorted.bam AddOrReplaceReadGroups.jar
	java $(javaArgs) -jar AddOrReplaceReadGroups.jar I=$< O=$@ LB=1 PL=illumina PU=whatever SM=whatever
	samtools index $@

marked.validated.aln.sorted.bam: validated.aln.sorted.bam MarkDuplicates.jar
	java $(javaArgs) -Djava.io.tmpdir=/tmp -jar MarkDuplicates.jar \
    INPUT=$< OUTPUT=$@ METRICS_FILE=metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

Chr4.dict: $(ref) CreateSequenceDictionary.jar
	java $(javaArgs) -jar CreateSequenceDictionary.jar R=$< O=$@

realigned.list: marked.validated.aln.sorted.bam $(ref) GenomeAnalysisTK.jar Chr4.dict
	java $(javaArgs) -jar GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R $(ref) -o $@ -I $<

realigned.bam: marked.validated.aln.sorted.bam realigned.list 
	java $(javaArgs) -Djava.io.tmpdir=/tmp  -jar GenomeAnalysisTK.jar -T IndelRealigner \
	-R $(ref) -o $@ -I $< -targetIntervals realigned.list

GenomeAnalysisTK.jar:
	ln -s /rhome/zwu/software/GATK/GenomeAnalysisTK.jar .
AddOrReplaceReadGroups.jar:
	ln -s /opt/picard/1.81/AddOrReplaceReadGroups.jar .
MarkDuplicates.jar:
	ln -s /opt/picard/1.81/MarkDuplicates.jar .
CreateSequenceDictionary.jar:
	ln -s /opt/picard/1.81/CreateSequenceDictionary.jar .

sam.raw.1.vcf gatk.raw.1.vcf: sam.raw.vcf gatk.raw.vcf 
	./generateFilteredFiles.py

plot.pdf venn_4.pdf: sam.raw.vcf gatk.raw.vcf varianttools.raw.vcf
	./overlap.py

venn.pdf: venn_4.pdf venn_1.pdf venn_2.pdf venn_3.pdf
	pdfjam --nup '2x2' --outfile $@ venn_1.pdf - venn_2.pdf - \
	venn_3.pdf - venn_4.pdf - 
plot.png:plot.pdf
	convert $< $@
venn.png:venn.pdf
	convert $< $@
