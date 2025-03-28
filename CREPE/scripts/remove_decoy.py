import pandas as pd

data = pd.read_csv('Clinvar1_offtarget_depth.txt',delimiter='\t',names=['chrom','pos','depth'])

chrom_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

data_clean = data[data['chrom'].isin(chrom_list)]

data_clean.to_csv('Clinvar1_offtarget_depth_nodecoy.txt',sep='\t',index=False)
