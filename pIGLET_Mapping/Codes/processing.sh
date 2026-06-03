#!/bin/bash

#environment needs samtools and seqtk


#fastq dir
fastq_dir="fastq"

#reference genome/full bam
reference_genome="/reference_genomes/hg38.fa"
full_bam="Full_Genome.bam"
sorted_full_bam="Full_Genome_sorted.bam"
full_fastq="Full_Genome.fastq"

#reference transgene/transgene bam
transgene_bam="Donor.bam"
sorted_transgene_bam="Donor_sorted.bam"
transgene_fasta="DonorSequence.fa"
transgene_mpileup="DonorSequence_sorted_mpileup_filter.txt"

#transgene read ids
begin_end_mpileup="begin_end_mpileup.txt"
qname_begin_end="qname_begin_end.txt"
qnames_newline="qnames_newline.txt"
qnames_uniq="qnames_uniq.txt"

#transgene reads aligned
transgene_ends_mapped="begin_end_qname_all_reads.fq"
transgene_aligned_bam="qname_begin_end_all_reads.bam"
transgene_aligned_bam_sorted="qname_begin_end_sorted_all_reads.bam"
transgene_aligned_sorted_sam="qname_begin_end_sorted_all_reads.sam"

#align the fastq reads and create a bam
dorado-1.2.0-linux-x64/bin/dorado aligner "$reference_genome" "$fastq_dir" --threads 32 > "$full_bam"

dorado-1.2.0-linux-x64/bin/dorado aligner "$transgene_fasta" "$fastq_dir" --threads 32 > "$transgene_bam"

cat "$fastq_dir"/*.fastq > "$full_fastq"

#sort and index bam file
#Program: samtools (Tools for alignments in the SAM format)

samtools sort -@ 32 "$full_bam" > "$sorted_full_bam"

samtools index  "$sorted_full_bam"


#sort and index transgene bam
samtools sort -@ 32 "$transgene_bam" > "$sorted_transgene_bam"

samtools index "$sorted_transgene_bam"

samtools faidx "$transgene_fasta"

#find coverage for every read
samtools-1.17/samtools mpileup -a -r Transgene_Fasta -f "$transgene_fasta" "$sorted_transgene_bam" -Q 20 -q 20 --output-QNAME > "$transgene_mpileup"


#extract sequences that align with the very beginning and end of the transgene
head -n 25 "$transgene_mpileup" > "$begin_end_mpileup"
tail -n 25 "$transgene_mpileup" >> "$begin_end_mpileup"


#extract the read names
awk -F'\t' '{print $7}' "$begin_end_mpileup" > "$qname_begin_end"
cat "$qname_begin_end" | tr "," "\n" > "$qnames_newline"
sort "$qnames_newline" | uniq > "$qnames_uniq"

#grab all sequences from the qnames uniq that are in the fastq
#cat together all fastq files created from the dorado basecaller step
#ends sequences that are mapping to the very ends of the transgene will also be mapped to the reference genome. Will be able to find exact location of transgene =
#Usage:   seqtk <command> <arguments>
#Version: 1.2-r94
#cat together all small fastqs into a larger fastq
seqtk subseq "$full_fastq" "$qnames_uniq" > "$transgene_ends_mapped"

/dorado-1.2.0-linux-x64/bin/dorado aligner "$reference_genome" "$transgene_ends_mapped" --threads 32 > "$transgene_aligned_bam"

samtools sort "$transgene_aligned_bam" > "$transgene_aligned_bam_sorted"
samtools index "$transgene_aligned_bam_sorted"
samtools view "$transgene_aligned_bam_sorted" > "$transgene_aligned_sorted_sam"


#look at transgene_aligned_sorted_sam output using cat transgene_aligned_sorted_sam.sam | less -S
#regions where many reads align to are where the transgene has been incorporated
