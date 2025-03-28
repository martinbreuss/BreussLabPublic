# The scripts in this directory were used to process and analyze targeted amplicon sequencing (TAS) data generated for primers created by CREPE. 

 tas_processing.py contains the pipeline to process TAS sequencing data from FASTQ to base recalibrated BAM file. pileup_maker.py creates pileup files for input target sites. tas_analysis.py uses the pileup files as input and calculates read depth, MAF, and the 95% CI for each input target site.

# Splitting On-target and Background Enrichment

After processing the FASTQ into a base recalibrated BAM file, this command was used to split the on-target coverage from the background enrichment, writing the on-target coverage into its own BAM file:

     samtools view -h -b -L NoOff.bed -o NoOff_target_sites.bam -@ 4 ClinVar1-96_S41_L003.recal.bam

This command was used to write the background data into a BAM file:

    samtools view -h -b -L NoOff.bed -U NoOff_background.bam -@ 4 ClinVar1-96_S41_L003.recal.bam


