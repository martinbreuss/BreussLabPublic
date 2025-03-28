import pandas as pd
import subprocess


data = pd.read_csv('libraries.csv', dtype=str)
#data = pd.read_csv('random96.csv')
#data = pd.read_csv('spotcheck.csv')

snp_loc = pd.read_csv('Clinvar1-96.csv', dtype=str)
#snp_loc = pd.read_csv('clinvar2_96_random.csv')
#snp_loc = pd.read_csv('clinvar2_spotcheck.csv')
#data = pd.read_csv('jareds.csv')

df_tas = pd.DataFrame(data)

TAS_list = df_tas.values.tolist()

#CHROM=['chr21', 'chr11', 'chr11', 'chr5', 'chr22', 'chr11', 'chr15', 'chr2', 'chr20', 'chr5', 'chr10', 'chr9', 'chr12', 'chr2']
#POS = ['45258259', '1172511', '119475487', '141441861', '20602713', '180318', '42528253', '100107183', '5066284', '139051170', '49739461', '35682258', '11120964', '31217463']
#ref_gt = ["T", "A", "T", "C", "T", "C", "C", "C", "G", "G", "C", "G", "G", "G", "G"]
#alt_gt = ["C", "C", "C", "T", "C", "T", "T", "G", "A", "A", "A", "A", "C", "C", "A"]

#CHROM=['chr6', 'chr16', 'chr19', 'chr3', 'chr4', 'chr10', 'chr4', 'chr2', 'chr17']
#POS = ['29268945', '79595584', '54055546', '179021806', '177961089', '124446608', '47406821', '170815568', '1727101']
#ref_gt = ["G", "A", "C", "C", "T", "G", "C", "C", "T"]
#alt_gt = ["A", "G", "G", "G", "C", "A", "T", "G", "C"]



for i in range(len(data['tas'])):
    for j in range(len(snp_loc['CHROM'])):
        subprocess.run(["samtools", "mpileup", "-r", 'chr'+str(snp_loc['CHROM'][j])+":"+str(snp_loc['POS'][j])+"-"+str(snp_loc['POS'][j]), "-f" ,"/pl/active/Breuss-Lab/reference_genomes/hg38.fa", "-Q", "15", "-q0", "-AB", "-d30000000", ""+str(data['tas'][i])+".recal.bam", "-o" ,"pileup/"+str(data['tas'][i])+"_"+str(snp_loc['CHROM'][j])+"_"+str(snp_loc['POS'][j])+".txt"])
