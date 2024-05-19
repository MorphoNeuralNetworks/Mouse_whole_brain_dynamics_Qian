import csv,os,json
import numpy as np
import nrrd
from scipy import stats
import seaborn as sns
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

with open('./Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

allen_exp=np.load("./Connectivity/Allen_connectivity_10.npy",allow_pickle=True)
allen_exp=allen_exp[0:len(TVMB_Re),:]
for i in range(0,len(TVMB_Re)): allen_exp[i,i]=0 
# allen_exp=np.log10(allen_exp+1)
# allen_exp[np.where(allen_exp<1e-5)]=0
allen_exp=allen_exp/np.amax(allen_exp)

connection=np.load("./Connectivity/SingleCell_connectivity_10.npy",allow_pickle=True)
connection=connection[0:len(TVMB_Re),:]
for i in range(0,len(TVMB_Re)): connection[i,i]=0 
# connection=np.log10(connection+1)
# connection[np.where(connection<1e-5)]=0
connection=connection/np.amax(connection)

area_color=np.load(r'D:\QPH\Data\Other_Infomation\color_network.npy', allow_pickle=True).item()
xtick=[x for x in TVMB_Re]+[x for x in TVMB_Re]
ytick=[x for x in TVMB_Re]+[x for x in TVMB_Re]
col_color=[area_color[x] for x in TVMB_Re]+[area_color[x] for x in TVMB_Re]
row_color=[area_color[x] for x in TVMB_Re]

import seaborn as sns
plt.close("all")
plot=sns.clustermap(data=allen_exp,figsize=(10,6),row_cluster=False,col_cluster=False,
               xticklabels=xtick,yticklabels=ytick,cmap=plt.get_cmap('hot_r'),#OrRd
               row_colors=row_color,col_colors=col_color,square=True,
               cbar_kws={'shrink':0.4,'label':"Connection strength"})
plt.tight_layout()
plt.title("Mesoscale from Allen")
plt.tight_layout()
# plt.savefig('Mesoscale from Allen.pdf', dpi=300)

plot=sns.clustermap(data=connection,figsize=(10,6),row_cluster=False,col_cluster=False,
               xticklabels=xtick,yticklabels=ytick,cmap=plt.get_cmap('hot_r'),#OrRd
               row_colors=row_color,col_colors=col_color,square=True,
               cbar_kws={'shrink':0.4,'label':"Connection strength"})
plt.tight_layout()
plt.title("Mesoscale from single neuron")
plt.tight_layout()
# plt.savefig('Mesoscale from single neuron.pdf', dpi=300)

l1=[]
l2=[]
for i in range(0,len(TVMB_Re)):
    for j in range(0,len(TVMB_Re)*2):
        if connection[i,j]!=0 and allen_exp[i,j]!=0:
            l2.append(connection[i,j])
            l1.append(allen_exp[i,j])

plt.close("all")
fig,ax1=plt.subplots(figsize=(4,3))
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

sns.kdeplot(l1,color='#1E90FF')
sns.kdeplot(l2,color='#DC143C')
ax1.set_xlabel("Normalized Connection Strength",size=14,fontproperties='Calibri',weight='bold')
ax1.set_ylabel("Density",size=14,fontproperties='Calibri',weight='bold')
x1_label = ax1.get_xticklabels() 
[x1_label_temp.set_fontname('Calibri') for x1_label_temp in x1_label]
y1_label = ax1.get_yticklabels() 
[y1_label_temp.set_fontname('Calibri') for y1_label_temp in y1_label]

ax2 = ax1.twinx()
sns.histplot(l1,bins=50,label="Allen",color='#1E90FF',ax=ax2,alpha=0.3,linewidth=0)
sns.histplot(l2,bins=50,label="single neuron",color='#DC143C',ax=ax2,alpha=0.3,linewidth=0)
ax2.set_ylabel("Count",fontsize=14,fontproperties='Calibri',weight='bold')
x1_label = ax2.get_xticklabels() 
[x1_label_temp.set_fontname('Calibri') for x1_label_temp in x1_label]
y1_label = ax2.get_yticklabels() 
[y1_label_temp.set_fontname('Calibri') for y1_label_temp in y1_label]
# plt.legend(prop={'family':'Calibri','weight':'bold','size':22},frameon=False,loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
plt.tight_layout()
# plt.savefig('Density of connection legend 5.pdf', dpi=300)

from scipy.stats import ks_2samp
import scipy.stats as stats
# s,p=stats.levene(l1, l2)
# if p>0.05:
#     print(stats.ttest_ind(l1, l2,equal_var=True))
# else:
#     print(stats.ttest_ind(l1, l2,equal_var=False))
print(len(l1))
print(ks_2samp(l1, l2))
print(stats.spearmanr(l1, l2))
# print(stats.pearsonr(l1, l2))
