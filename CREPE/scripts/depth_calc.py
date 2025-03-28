import pandas as pd
import numpy as np

data = pd.read_csv('Clinvar1_target_depth.txt',delimiter='\t',names=['chrom','pos','depth'])

data_nozero = data[data['depth'] != 0]

#Clinvar1_target_depth.txt
#Clinvar1_offtarget_depth.txt

mean_depth = np.mean(data_nozero['depth'])

max_value = data_nozero.max()

print(len(data_nozero['depth']))
print(max_value)
print(mean_depth)

data_high = data_nozero[data_nozero['depth'] > 538]

mean_high = np.mean(data_high['depth'])

max_high = data_high.max()


print(len(data_high['depth']))
print(mean_high)
print(max_high)


