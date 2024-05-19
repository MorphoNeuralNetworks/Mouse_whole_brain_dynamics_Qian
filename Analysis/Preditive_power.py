import numpy as np
import time,shutil,csv,random
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

color_dict=np.load(r"..\Data\color_network.npy",allow_pickle=True).item()
c_list=[color_dict[i] for i in TVMB_Re]

fcd_exp_all=np.load("FCD_exp.npy",allow_pickle=True)
kk=[len(x) for x in fcd_exp_all]

count=len(fcd_exp_all)
ii=1
plt.close("all")
# plt.figure(figsize=(16,4))
single_par,allen_par=np.load("best_par.npy",allow_pickle=True)
for i in range(len(allen_par)):
    if allen_par[i,0]<0.084:
        allen_par[i,0]=0.084
'''
# 计算FCD的Predictive Power
pp_list=[]
fcd_exp_all=np.load("FCD_exp.npy",allow_pickle=True)
for exp_id in range(0,len(fcd_exp_all)):
    
    fcd_exp=fcd_exp_all[exp_id]
    triangular = np.triu_indices(len(fcd_exp), 1)
    tr_exp=fcd_exp[triangular]
    x_exp, y_exp = sns.kdeplot(tr_exp).get_lines()[0].get_data()
    plt.close()
    
    padding_s,padding_e=60,660
        
    filepath="./FCD_single_detail/"
    coupling_v,noise_v=single_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd_single=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    triangular = np.triu_indices(len(fcd_single), 1)
    fcd_single_up=fcd_single[triangular]
    x_single, y_single = sns.kdeplot(fcd_single_up).get_lines()[0].get_data()
    plt.close()
   
    filepath="./FCD_allen_detail/"
    coupling_v,noise_v=allen_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd_allen=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    triangular = np.triu_indices(len(fcd_allen), 1)
    fcd_allen_up=fcd_allen[triangular]
    x_allen, y_allen = sns.kdeplot(fcd_allen_up).get_lines()[0].get_data()
    plt.close()
    # plt.legend(["Exp","Single","Allen"])
    
    pp_list.append([np.sum(np.absolute(y_single-y_exp))+np.sum(np.absolute(x_single-x_exp)), 
                    np.sum(np.absolute(y_allen-y_exp))+np.sum(np.absolute(x_allen-x_exp))])

# np.save("pp_list.npy",pp_list)
'''

# 平均BOLD的预测结果
pp_list=[]
bold_data_exp_all=np.load("bold_data_exp.npy",allow_pickle=True)
for exp_id in range(0,len(bold_data_exp_all)):
    bold_data_exp=bold_data_exp_all[exp_id]
    for i in range(len(bold_data_exp)):
        bold_data_exp[i,:]=(bold_data_exp[i,:]-np.min(bold_data_exp[i,:]))/(np.max(bold_data_exp[i,:])-np.min(bold_data_exp[i,:]))
        # bold_data_exp[i,:]=(bold_data_exp[i,:]-np.mean(bold_data_exp[i,:]))/(np.std(bold_data_exp[i,:]))
    padding_s,padding_e=60,660
        
    filepath="./BOLD_single_detail/"
    coupling_v,noise_v=single_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    bold_data_single=np.load(filepath+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data_single=bold_data_single[:,0,:,0].T
    bold_data_single=bold_data_single[:,padding_s:padding_e]
    bold_data_single=(bold_data_single-np.min(bold_data_single))/(np.max(bold_data_single)-np.min(bold_data_single))
   
    filepath="./BOLD_allen_detail/"
    coupling_v,noise_v=allen_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    bold_data_allen=np.load(filepath+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data_allen=bold_data_allen[:,0,:,0].T
    bold_data_allen=bold_data_allen[:,padding_s:padding_e]
    bold_data_allen=(bold_data_allen-np.min(bold_data_allen))/(np.max(bold_data_allen)-np.min(bold_data_allen))
    
    #皮层区域
    ave_single=np.mean(bold_data_single[0:14,:],1)
    ave_allen=np.mean(bold_data_allen[0:14,:],1)
    ave_exp=np.mean(bold_data_exp[0:14,:],1)
    #丘脑区域
    # ave_single=np.mean(bold_data_single[19:28,:],1)
    # ave_allen=np.mean(bold_data_allen[19:28,:],1)
    # ave_exp=np.mean(bold_data_exp[19:28,:],1)
    pp_list.append((np.sum(np.absolute(ave_single-ave_exp)),
                    np.sum(np.absolute(ave_allen-ave_exp))))

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# pp_list=np.load("pp_list.npy",allow_pickle=True)
# for i in [17,16,12,10,1,0]:
#     pp_list=np.delete(pp_list,i,0)
plt.figure(figsize=(16,2))
grid=plt.GridSpec(1,8)
plt.subplot(grid[0,0:6])
width=0.4
x_t=np.array([x for x in range(1,len(pp_list)+1)])
pp_list=np.array(pp_list)
plt.bar(x=range(1,len(pp_list)+1),height=pp_list[:,0],width=0.4,color='#FFB6C1')
plt.bar(x=x_t+width,height=pp_list[:,1],width=0.4,color='#87CEFA')
plt.xticks([x for x in range(1,len(pp_list)+1)])
plt.xlabel("Exp id", fontsize=14)
plt.ylabel("Loss", fontsize=14)
plt.legend(["Single","Allen"],prop={'family':'Calibri','weight':'bold','size':18},frameon=False)

plt.subplot(grid[0,6:8])
bplot=plt.boxplot(pp_list,vert=True,patch_artist=True)

colors = ['#FFB6C1', '#87CEFA', 'lightgreen']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

plt.xticks([1,2],["Single","Allen"], fontsize=14)
plt.ylabel("Loss", fontsize=14)
# plt.title("Loss of FCD distirbution", fontsize=18)
plt.tight_layout()
# plt.savefig("./Loss bar th.pdf",dpi=300)

from scipy.stats import ks_2samp
import scipy.stats as stats
def myttest(A,B):
    import scipy.stats as stats
    s,p=stats.levene(A,B)
    if p>0.05:
        print(stats.ttest_ind(A,B,equal_var=True))
    else:
        print(stats.ttest_ind(A,B,equal_var=False))

print(myttest(pp_list[:,0],pp_list[:,1]))
print(ks_2samp(pp_list[:,0],pp_list[:,1]))
from scipy.stats import shapiro

stat, p = shapiro(pp_list[:,0])
print((stat, p))
if p > 0.05:
    print('不能拒绝原假设，样本数据服从正态分布')
else:
    print('不服从正态分布')
stat, p = shapiro(pp_list[:,1])
print((stat, p))
if p > 0.05:
    print('不能拒绝原假设，样本数据服从正态分布')
else:
    print('不服从正态分布')



