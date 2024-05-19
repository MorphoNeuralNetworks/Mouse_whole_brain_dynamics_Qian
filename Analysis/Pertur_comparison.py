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

connection_o=np.load("./Connectivity/SingleCell_connectivity_10.npy",allow_pickle=True)
connection_o=connection_o[0:len(TVMB_Re),:]
connection_o=np.log10(connection_o+1)
connection=connection_o/np.amax(connection_o)
connection=connection_o

scale=np.load("./Connectivity/Perturbation/scale_0.5_prune_0_delete_0_all_SingleCell_connectivity_10.npy",allow_pickle=True)
scale=scale[0:len(TVMB_Re),:]
scale=np.log10(scale+1)
scale=scale/np.amax(connection_o)

prune=np.load("./Connectivity/Perturbation/scale_0_prune_0.5_delete_0_all_SingleCell_connectivity_10.npy",allow_pickle=True)
prune=prune[0:len(TVMB_Re),:]
prune=np.log10(prune+1)
prune=prune/np.amax(connection_o)
delete=np.load("./Connectivity/Perturbation/scale_0_prune_0_delete_0.5_all_SingleCell_connectivity_10.npy",allow_pickle=True)
delete=delete[0:len(TVMB_Re),:]
delete=np.log10(delete+1)
delete=delete/np.amax(connection_o)


l1,l2,l3,l4=[],[],[],[]
for i in range(len(TVMB_Re)):
    for j in range(len(TVMB_Re)*2):
        if connection[i,j]!=0:
            l1.append(connection[i,j])
        if scale[i,j]!=0:
            l2.append(scale[i,j])
        if prune[i,j]!=0:
            l3.append(prune[i,j])
        if delete[i,j]!=0:
            l4.append(delete[i,j])
from scipy.stats import ks_2samp
import scipy.stats as stats
print(ks_2samp(l1, l2))
# print(stats.spearmanr(l1, l2))

print(ks_2samp(l1, l3))
# print(stats.spearmanr(l1, l3))

print(ks_2samp(l1, l4))
# print(stats.spearmanr(l1, l4))

plt.close("all")
fig,ax1=plt.subplots(figsize=(4,3))
ax1.spines['top'].set_visible(False)
sns.kdeplot(l1,color='#1f77b4',label="Non-Perturbed",linewidth=3)
sns.kdeplot(l2,color='#d62728',label="Scale 0.5",linewidth=3)
sns.kdeplot(l3,color='#2ca02c',label="Prune 0.5",linewidth=3)
sns.kdeplot(l4,color='#9467bd',label="Delete 0.5",linewidth=3)
ax1.set_xlabel("Normalized Connection Strength",size=14,fontproperties='Calibri',weight='bold')
ax1.set_ylabel("Density",size=14,fontproperties='Calibri',weight='bold')
x1_label = ax1.get_xticklabels() 
[x1_label_temp.set_fontname('Calibri') for x1_label_temp in x1_label]
y1_label = ax1.get_yticklabels() 
[y1_label_temp.set_fontname('Calibri') for y1_label_temp in y1_label]
# plt.legend()
ax2 = ax1.twinx()
ax2.spines['top'].set_visible(False)
sns.histplot(l1,bins=50,label="Non-Perturbed",color='#1f77b4',ax=ax2,alpha=1,linewidth=0)
sns.histplot(l2,bins=50,label="Scale 0.5",color='#d62728',ax=ax2,alpha=0.8,linewidth=0)
sns.histplot(l3,bins=50,label="Prune 0.5",color='#2ca02c',ax=ax2,alpha=0.8,linewidth=0)
sns.histplot(l4,bins=50,label="Delete 0.5",color='#9467bd',ax=ax2,alpha=0.8,linewidth=0)
ax2.set_ylabel("Count",fontsize=14,fontproperties='Calibri',weight='bold')
x1_label = ax2.get_xticklabels() 
[x1_label_temp.set_fontname('Calibri') for x1_label_temp in x1_label]
y1_label = ax2.get_yticklabels() 
[y1_label_temp.set_fontname('Calibri') for y1_label_temp in y1_label]

# plt.legend(prop={'family':'Calibri','weight':'bold','size':22},frameon=False,loc=2, bbox_to_anchor=(1.05,1.0),borderaxespad = 0.)
plt.tight_layout()
# plt.savefig('Density of perturbation connection.pdf', dpi=300)
'''
'''
'''
area_color=np.load(r'D:\QPH\Data\Other_Infomation\color_network.npy', allow_pickle=True).item()
xtick=[x for x in TVMB_Re]+[x for x in TVMB_Re]
ytick=[x for x in TVMB_Re]+[x for x in TVMB_Re]
col_color=[area_color[x] for x in TVMB_Re]+[area_color[x] for x in TVMB_Re]
row_color=[area_color[x] for x in TVMB_Re]

plt.close("all")
plot=sns.clustermap(data=connection,figsize=(10,6),row_cluster=False,col_cluster=False,
               xticklabels=xtick,yticklabels=ytick,cmap=plt.get_cmap('hot_r'),#OrRd
               row_colors=row_color,col_colors=col_color,square=True,vmin=0,vmax=1,
               cbar_kws={'shrink':0.4,'label':"Connection strength"})
plt.tight_layout()
plt.title("Non-Perturbed")
plt.tight_layout()
plt.savefig('Non-Perturbed.pdf', dpi=300)

plot=sns.clustermap(data=delete,figsize=(10,6),row_cluster=False,col_cluster=False,
               xticklabels=xtick,yticklabels=ytick,cmap=plt.get_cmap('hot_r'),#OrRd
               row_colors=row_color,col_colors=col_color,square=True,vmin=0,vmax=1,
               cbar_kws={'shrink':0.4,'label':"Connection strength"})
plt.tight_layout()
plt.title("Delete 0.5")
plt.tight_layout()
plt.savefig('Delete 0.5.pdf', dpi=300)
'''