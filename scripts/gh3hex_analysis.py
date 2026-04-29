#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 17:55:26 2025

@author: Sebastian
"""




import glob
import pandas as pd 
import os 
import matplotlib.pyplot as plt  
import seaborn as sns

from scipy.stats import sem
import numpy as np
import pingouin as pg

import locale 

#### GH3 related mutants 
sam_size = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/gh3hex_sam_size.xlsx')
sam_size["length"] = np.sqrt(sam_size['area']/3.14)*2


plt.figure(figsize=(1,3))
ax = sns.stripplot(x="genotype", y='length', data = sam_size, 

                  order= ['Col0', 'gh3hex'], 
                  palette = 'inferno',  dodge = False,
                  linewidth=0, size = 10,  alpha  = 0.3, 
                   zorder=0)
ax = sns.pointplot(x="genotype", y='length', 
                  data = sam_size, 
                  # linecolor = 'black', 
                  order= ['Col0', 'gh3hex'], 
                   join=False,
                   capsize=.2,  palette = ("black", "black"), scale=0.7,
                   errwidth = 1
                   )

ax.set_ylim(0,120)
plt. grid(False)
plt.setp(ax.artists, edgecolor = 'k', facecolor='black')
plt.setp(ax.lines, color='k')
sns.despine(offset=5, trim=False)
ax.legend(bbox_to_anchor=(2, 1.1))
plt.xticks(rotation=45, fontsize = 10)
ax.set_ylabel("SAM size (µm)" , fontsize =  12)
ax.set_xlabel("" , fontsize =  12)
plt.xticks([0, 1],  ['Col-0',  'gh3hex'], fontsize = 12,rotation = 20)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/gh3hex_sam_size.pdf',  bbox_inches='tight')
plt.show()


# temp1 = amidrc[amidrc["nitrate"] == nit]
result = pg.pairwise_tukey(data = sam_size, dv='length', between=['genotype']).round(10)
print(result)


### DIVERGE ANGLE FOR gh3hex mutant lines 
angle = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/diverge_angle.xlsx')
angle['plant'] = angle['plant'] + '_' +  angle['genotype']


def calculate_divergence_angle(df):
    angles = df['angle'].values
    divergence_angles = []
    
    for i in range(len(angles) - 1):
        diff = angles[i+1] - angles[i]
        if diff < 0:
            diff += 360
        divergence_angles.append(diff)
    
    divergence_angles.append(None)  # The last row has no next row to compare with
    df['angle_final'] = divergence_angles
    return df



angle_final = pd.DataFrame()
for plant in angle['plant'].unique(): 
    temp1 = angle[angle['plant'] == plant]
    df_with_divergence = calculate_divergence_angle(temp1)
    angle_final = pd.concat([angle_final,df_with_divergence])
    
angle_final.to_excel('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/diverge_angle_final.xlsx')
angle_final = angle_final.dropna(subset=['angle_final'])

    
plt.figure(figsize=(4,2.5))
sns.histplot(data=angle_final, x='angle_final', kde=True, hue  = 'genotype',
             stat="count", palette = ('purple', 'black'), label="Probabilities")
plt.xlabel('Divergence angles', fontsize = 12)
plt.ylabel('Frequency of siliques', fontsize = 12)
plt.grid(False)
plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/gh3hex_phylotaxis.pdf',  bbox_inches='tight')
plt.show()



# ############################################# 
# angle = pd.read_excel('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/diverge_angle.xlsx')
# angle['replicate'] = angle.plant.str.split('p', expand = True)[0]
# angle = angle[angle['replicate'] != 'R1']

# def calculate_divergence_angle(df):
#     angles = df['angle'].values
#     divergence_angles = []
    
#     for i in range(len(angles) - 1):
#         diff = angles[i+1] - angles[i]
#         if diff < 0:
#             diff += 360
#         divergence_angles.append(diff)
    
#     divergence_angles.append(None)  # The last row has no next row to compare with
#     df['angle_final'] = divergence_angles
#     return df


# angle_final = pd.DataFrame()

# for plant in angle['plant'].unique(): 
#     temp1 = angle[angle['plant'] == plant]
#     df_with_divergence = calculate_divergence_angle(temp1)
#     angle_final = pd.concat([angle_final,df_with_divergence ])
    
       
# plt.figure(figsize=(4,2.5))
# sns.lineplot(data=angle_final, x='silique',
#              y = 'angle_final', 
#             hue  = 'genotype',
#              palette = ('purple', 'black'))
# plt.xlabel('Silique (from base to top)', fontsize = 10)
# plt.ylabel('Diverge angle', fontsize = 10)
# plt.grid(False)
# # plt.savefig('/Users/Sebastian/Documentos/SLCU_lab/results/scrna_seq/validation_experiments/gh3hex/gh3hex_phylotaxis.pdf',  bbox_inches='tight')
# plt.show()

from scipy import stats




pg.anderson(x =angle_final[angle_final['genotype'] =='col0']['angle_final'].tolist(),
            y = angle_final[angle_final['genotype'] =='gh3hex']['angle_final'].tolist())

pg.normality( data =angle_final , dv = 'angle_final', group = 'genotype')

stats.ks_2samp( angle_final[angle_final['genotype'] =='col0']['angle_final'].tolist(), 
               angle_final[angle_final['genotype'] =='gh3hex']['angle_final'].tolist(),
               )











