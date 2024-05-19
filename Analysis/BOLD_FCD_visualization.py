import numpy as np
import time,shutil,csv,random,os,shutil
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
from scipy import signal

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

color_dict=np.load(r"..\Data\color_network.npy",allow_pickle=True).item()
c_list=np.array([color_dict[i] for i in TVMB_Re])

fcd_exp_all=np.load("FCD_exp.npy",allow_pickle=True)
bold_data_exp_all=np.load("bold_data_exp.npy",allow_pickle=True)

count=len(fcd_exp_all)
ii=1

single_par,allen_par=np.load("best_par.npy",allow_pickle=True)
for i in range(len(allen_par)):
    if allen_par[i,0]<0.084:
        allen_par[i,0]=0.084
print(np.mean(single_par,0))
print(np.std(single_par,0))
print(np.mean(allen_par,0))
print(np.std(allen_par,0))


writepath="./Figures/"
if os.path.exists(writepath):
    shutil.rmtree(writepath)  
os.mkdir(writepath) 

plt.close("all")
count=0
loss=[]
for exp_id in range(6,7):#len(bold_data_exp_all)):

    bold_data_exp=bold_data_exp_all[exp_id]
    for i in range(len(bold_data_exp)):
        bold_data_exp[i,:]=(bold_data_exp[i,:]-np.min(bold_data_exp[i,:]))/(np.max(bold_data_exp[i,:])-np.min(bold_data_exp[i,:]))
    # bold_data_exp=(bold_data_exp-np.min(bold_data_exp))/(np.max(bold_data_exp)-np.min(bold_data_exp))
    fcd_exp=fcd_exp_all[exp_id]
    triangular = np.triu_indices(len(fcd_exp), 1)
    fcd_exp_up=fcd_exp[triangular]
    padding_s,padding_e=60,660
        
    filepath="./BOLD_single/"
    coupling_v,noise_v=single_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    bold_data_single=np.load(filepath+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data_single=bold_data_single[:,0,:,0].T
    bold_data_single=bold_data_single[:,padding_s:padding_e]
    
    # for i in range(len(bold_data_single)):
    #     bold_data_single[i,:]=(bold_data_single[i,:]-np.min(bold_data_single[i,:]))/(np.max(bold_data_single[i,:])-np.min(bold_data_single[i,:]))
    bold_data_single=(bold_data_single-np.min(bold_data_single))/(np.max(bold_data_single)-np.min(bold_data_single))
    filepath="./FCD_single/"
    coupling_v,noise_v=single_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd_single=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    triangular = np.triu_indices(len(fcd_single), 1)
    fcd_s_up=fcd_single[triangular]
   
    filepath="./BOLD_allen/"
    coupling_v,noise_v=allen_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    bold_data_allen=np.load(filepath+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data_allen=bold_data_allen[:,0,:,0].T
    bold_data_allen=bold_data_allen[:,padding_s:padding_e]
    
    # for i in range(len(bold_data_single)):
    #     bold_data_allen[i,:]=(bold_data_allen[i,:]-np.min(bold_data_allen[i,:]))/(np.max(bold_data_allen[i,:])-np.min(bold_data_allen[i,:]))
    bold_data_allen=(bold_data_allen-np.min(bold_data_allen))/(np.max(bold_data_allen)-np.min(bold_data_allen))
    filepath="./FCD_allen/"
    coupling_v,noise_v=allen_par[exp_id]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd_allen=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    triangular = np.triu_indices(len(fcd_allen), 1)
    fcd_a_up=fcd_allen[triangular]
    
    
    bold_time=[x for x in range(1,601)]
    select_id=[22,17,5,24,0]
    num=len(TVMB_Re)
    plt.close("all")
    plt.figure(figsize=(16,10))
    grid=plt.GridSpec(3,4)
    # BOLD exp
    plt.subplot(grid[0,0:2])
    for k in range(0,num):
        if k in select_id:
            plt.plot(bold_time,bold_data_exp[k,:],color=c_list[k],linewidth=1.5)
        else:
            plt.plot(bold_time,bold_data_exp[k,:],color=c_list[k],linewidth=0.5,alpha=0.3)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Exp', fontsize=14)
    # barplot
    plt.subplot(grid[0,2])
    bold_average_exp=np.mean(bold_data_exp,1)
    err_exp=np.std(bold_data_exp,1)
    # plt.bar(range(0,num),bold_average_exp[0:num],color=c_list,alpha=0.3,yerr=err_exp[0:num],error_kw={'ecolor':"0.1",'capsize':2,'alpha':0.3})
    # plt.bar(select_id,bold_average_exp[select_id],color=c_list[select_id],alpha=1,yerr=err_exp[select_id],error_kw={'ecolor':"0.1",'capsize':2})
    plt.bar(range(0,num),bold_average_exp[0:num],color=c_list,yerr=err_exp[0:num],error_kw={'ecolor':"0.1",'capsize':2})
    plt.xticks(range(0,num), TVMB_Re, fontsize=8,rotation=270)
    plt.ylabel('Average Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Exp', fontsize=14)
    #FCD
    plt.subplot(grid[0,3])
    cs=plt.imshow(fcd_exp, cmap='jet',vmin=0,vmax=1) 
    axcb=plt.colorbar(ticks=[0, 0.5, 1])
    cs.set_clim(0, 1.0)
    plt.xlabel(r'Time $t_j$ (s)', fontsize=12)
    plt.ylabel(r'Time $t_i$ (s)', fontsize=12)
    plt.title('FCD Exp', fontsize=14)
    
    # BOLD single
    plt.subplot(grid[1,0:2])
    for k in range(0,num):
        if k in select_id:
            plt.plot(bold_time,bold_data_single[k,:],color=c_list[k],linewidth=1.5)
        else:
            plt.plot(bold_time,bold_data_single[k,:],color=c_list[k],linewidth=0.5,alpha=0.3)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Single', fontsize=14)
    # barplot
    plt.subplot(grid[1,2])
    bold_average_single=np.mean(bold_data_single,1)
    err_single=np.std(bold_data_single,1)
    # plt.bar(range(0,num),bold_average_single[0:num],color=c_list,alpha=0.3,yerr=err_single[0:num],error_kw={'ecolor':"0.1",'capsize':2,'alpha':0.3})
    # plt.bar(select_id,bold_average_single[select_id],color=c_list[select_id],alpha=1,yerr=err_single[select_id],error_kw={'ecolor':"0.1",'capsize':2})
    plt.bar(range(0,num),bold_average_single[0:num],color=c_list,yerr=err_single[0:num],error_kw={'ecolor':"0.1",'capsize':2})
    plt.xticks(range(0,num), TVMB_Re, fontsize=8,rotation=270)
    plt.ylabel('Average Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Single', fontsize=14)
    #FCD
    plt.subplot(grid[1,3])
    cs=plt.imshow(fcd_single, cmap='jet',vmin=0,vmax=1) 
    axcb=plt.colorbar(ticks=[0, 0.5, 1])
    cs.set_clim(0, 1.0)
    plt.xlabel(r'Time $t_j$ (s)', fontsize=12)
    plt.ylabel(r'Time $t_i$ (s)', fontsize=12)
    plt.title('FCD Single', fontsize=14)
    
    # BOLD allen
    plt.subplot(grid[2,0:2])
    for k in range(0,num):
        if k in select_id:
            plt.plot(bold_time,bold_data_allen[k,:],color=c_list[k],linewidth=1.5)
        else:
            plt.plot(bold_time,bold_data_allen[k,:],color=c_list[k],linewidth=0.5,alpha=0.3)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Allen', fontsize=14)
    # barplot
    plt.subplot(grid[2,2])
    bold_average_allen=np.mean(bold_data_allen,1)
    err_allen=np.std(bold_data_allen,1)
    # plt.bar(range(0,num),bold_average_allen[0:num],color=c_list,alpha=0.3,yerr=err_allen[0:num],error_kw={'ecolor':"0.1",'capsize':2,'alpha':0.3})
    # plt.bar(select_id,bold_average_allen[select_id],color=c_list[select_id],alpha=1,yerr=err_allen[select_id],error_kw={'ecolor':"0.1",'capsize':2})
    plt.bar(range(0,num),bold_average_allen[0:num],color=c_list,yerr=err_allen[0:num],error_kw={'ecolor':"0.1",'capsize':2})
    plt.xticks(range(0,num), TVMB_Re, fontsize=8,rotation=270)
    plt.ylabel('Average Amplitude (%)', fontsize=12)
    plt.title('BOLD Right Allen', fontsize=14)
    #FCD
    plt.subplot(grid[2,3])
    cs=plt.imshow(fcd_allen, cmap='jet',vmin=0,vmax=1) 
    axcb=plt.colorbar(ticks=[0, 0.5, 1])
    cs.set_clim(0, 1.0)
    plt.xlabel(r'Time $t_j$ (s)', fontsize=12)
    plt.ylabel(r'Time $t_i$ (s)', fontsize=12)
    plt.title('FCD Alllen', fontsize=14)

    plt.tight_layout()
    plt.savefig("./Figures/"+str(exp_id)+".pdf",dpi=300)
    # plt.savefig("./Figures/"+str(exp_id)+".jpg",dpi=300)

