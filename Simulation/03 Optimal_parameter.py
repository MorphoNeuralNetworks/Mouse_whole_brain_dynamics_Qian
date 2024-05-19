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

# 寻找每个实验对应的最优参数组合
allen_par=[]
single_par=[]
for exp_id in range(0,len(fcd_exp_all)):
    print(exp_id)
    fcd_exp=fcd_exp_all[exp_id]
    triangular = np.triu_indices(len(fcd_exp), 1)
    fcd_exp_up=fcd_exp[triangular]
    k1=sns.kdeplot(fcd_exp_up) 
    x1, y1 = sns.kdeplot(fcd_exp_up).get_lines()[0].get_data()
    plt.close()
    
    loss_matrix=np.zeros((100,100))
    count=0
    filepath='./kdeplot_single/'
    for i in range(0,100):
        for j in range(0,100):
            count+=1
            (coupling_v,noise_v)=(i*0.001+0.2,j*0.000001+0.000001)
            file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
            [x2, y2]=np.load(filepath+"kdeplot"+file_end+".npy",allow_pickle=True)
            loss=np.sum(np.absolute(y2-y1))+np.sum(np.absolute(x2-x1)) #绝对值
            loss_matrix[i,j]=loss
    print(np.where(loss_matrix==np.min(loss_matrix)))
    [x_t,y_t]=np.where(loss_matrix==np.min(loss_matrix))
    single_par.append((x_t[0]*0.001+0.2,y_t[0]*0.000001+0.000001))
    
    loss_matrix=np.zeros((100,200))
    count=0
    filepath='./kdeplot_allen/'
    for i in range(0,100):
        for j in range(0,200):
            count+=1
            (coupling_v,noise_v)=(i*0.001+0.05,j*0.000001+0.000001)
            file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
            [x2, y2]=np.load(filepath+"kdeplot"+file_end+".npy",allow_pickle=True)
            
            loss=np.sum(np.absolute(y2-y1))+np.sum(np.absolute(x2-x1)) #绝对值
            loss_matrix[i,j]=loss
    print(np.where(loss_matrix==np.min(loss_matrix)))
    [x_t,y_t]=np.where(loss_matrix==np.min(loss_matrix))
    allen_par.append((x_t[0]*0.001+0.05,y_t[0]*0.000001+0.000001))
np.save("best_par.npy",[single_par,allen_par])


