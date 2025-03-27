import subprocess
import pandas as pd

#filepath = '/scratch/alpine/jonpitsch@xsede.org/tas/sperm3_tas/'

filepath = '/scratch/alpine/jonpitsch@xsede.org/tas/clinvar2/' 

data = pd.read_csv('/scratch/alpine/jonpitsch@xsede.org/tas/clinvar2/libraries.csv')

df = pd.DataFrame(data)

TAS_list = df.values.tolist()

for i in range(len(TAS_list)):
    subprocess.run(('bwa mem -t 20 /pl/active/Breuss-Lab/reference_genomes/hg38.fa '+filepath+df['tas'][i]+'_R1_001.fastq.gz '+filepath+df['tas'][i]+'_R2_001.fastq.gz| samtools view -@ 20 -Sb - > '+filepath+df['tas'][i]+'.bam'), shell=True)

for i in range(len(TAS_list)):
    subprocess.run(('samtools sort -@ 20 -o '+filepath+df['tas'][i]+'.sorted.bam '+filepath+df['tas'][i]+'.bam'), shell=True)

for i in range(len(TAS_list)):
    subprocess.run(('samtools index -@ 20 '+filepath+df['tas'][i]+'.sorted.bam'), shell=True)

for i in range(len(TAS_list)):
    subprocess.run(('picard -Xmx60G FixMateInformation I='+filepath+df['tas'][i]+'.sorted.bam O='+filepath+df['tas'][i]+'.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR='+filepath+' MAX_RECORDS_IN_RAM=50000'), shell=True)

for i in range(len(TAS_list)):
    subprocess.run(('picard -Xmx60G AddOrReplaceReadGroups I='+filepath+df['tas'][i]+'.fixed.bam O='+filepath+df['tas'][i]+'.add.bam LB=lib1 PL=ILLUMINA PU=unit1 SM='+df['tas'][i]+' CREATE_INDEX=true TMP_DIR='+filepath+' MAX_RECORDS_IN_RAM=50000'), shell=True)


for i in range(len(TAS_list)):
    subprocess.run(('/projects/$USER/DeepMosaic/gatk-4.2.6.1/gatk BaseRecalibrator -R /pl/active/Breuss-Lab/reference_genomes/hg38.fa -I '+filepath+df['tas'][i]+'.add.bam --known-sites /pl/active/Breuss-Lab/gatk_resources/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /pl/active/Breuss-Lab/gatk_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf --known-sites /pl/active/Breuss-Lab/gatk_resources/Homo_sapiens_assembly38.known_indels.vcf -O '+filepath+df['tas'][i]+'.recal.table'), shell=True)


for i in range(len(TAS_list)):
    subprocess.run(('/projects/$USER/DeepMosaic/gatk-4.2.6.1/gatk ApplyBQSR -R /pl/active/Breuss-Lab/reference_genomes/hg38.fa -I '+filepath+df['tas'][i]+'.add.bam --bqsr-recal-file '+filepath+df['tas'][i]+'.recal.table -O '+filepath+df['tas'][i]+'.recal.bam'), shell=True)
