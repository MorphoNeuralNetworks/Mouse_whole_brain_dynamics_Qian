import numpy as np
import time,shutil,csv,random
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

color_dict=np.load(r"..\Data\color_network.npy",allow_pickle=True).item()
c_list=[color_dict[i] for i in TVMB_Re]

count=0
loss_matrix=np.zeros((100,200))

plt.close("all")
fcd_exp_all=np.load("FCD_exp.npy",allow_pickle=True)
kk=[len(x) for x in fcd_exp_all]

for exp_id in range(0,20):
    fcd_exp=fcd_exp_all[exp_id]
    triangular = np.triu_indices(len(fcd_exp), 1)
    fcd_exp_up=fcd_exp[triangular]
    
    k1=sns.kdeplot(fcd_exp_up) 
    x1, y1 = sns.kdeplot(fcd_exp_up).get_lines()[0].get_data()
    # plt.close()
    count=0
    filepath='./kdeplot_single/'

    for i in range(0,100):
        for j in range(0,100):
            count+=1
            (coupling_v,noise_v)=(i*0.001+0.2,j*0.000001+0.000001)      
            file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
            [x2, y2]=np.load(filepath+"kdeplot"+file_end+".npy",allow_pickle=True)
            # loss=np.sum((y2-y1)**2) #平方和
            # loss_matrix[i,j]=loss
            
            # loss=np.sum(np.absolute(y2-y1)) #绝对值
            # loss_matrix[i,j]=loss
            
            # loss=np.sum((y2-y1)**2)+np.sum((x2-x1)**2) #欧氏距离
            # loss_matrix_t1[i,j]=loss
            
            loss=np.sum(np.absolute(y2-y1))+np.sum(np.absolute(x2-x1)) #曼哈顿距离
            loss_matrix[i,j]=loss
    plt.close("all")
    plt.figure(figsize=(6,4))
    cs=plt.imshow(loss_matrix, cmap='jet') 
    axcb=plt.colorbar()
    axcb.set_label('Loss', fontsize=14)
    plt.yticks([0,50,99],['0.2','0.25','0.3'], fontsize=12)
    plt.xticks([0,50,99],['1e-6','5e-5', '1e-4'], fontsize=12)

    plt.xlabel(r'Coupling', fontsize=14)
    plt.ylabel(r'Noise', fontsize=14)
    plt.tight_layout()
    # plt.savefig("./Figures/Loss_"+str(exp_id)+".jpg",dpi=300)
    # print(np.where(loss_matrix==np.max(loss_matrix)))