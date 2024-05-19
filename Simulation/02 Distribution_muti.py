import numpy as np
import time,shutil
import os
import seaborn as sns
import matplotlib.pyplot as plt
name="single"
def Get_kdeplot(par):
    writepath="./kdeplot_"+name+"/"
    filepath='./FCD_'+name+'/'
    coupling_v,noise_v=par[0],par[1]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    triangular = np.triu_indices(len(fcd), 1)
    fcd_up=fcd[triangular]
    
    x, y = sns.kdeplot(fcd_up).get_lines()[0].get_data()
    plt.close()
    np.save(writepath+"kdeplot"+file_end+".npy",[x,y])
    
def run__pool():  # main process

    from multiprocessing import Pool
    cpu_worker_num = 36 #设定并行的数量
    writepath="./kdeplot_"+name+"/"
    if os.path.exists(writepath):
        shutil.rmtree(writepath)  
    os.mkdir(writepath) 

    par_list=[]
    # set the parameters range
    for i in range(0,100):
        for j in range(0,100):
            par_list.append((i*0.001+0.2,j*0.000001+0.000001))
 
    import time
    time_start = time.time()  # 记录开始时间
    with Pool(cpu_worker_num) as p:
        p.map(Get_kdeplot, par_list)
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('Get Dis time: '+str(time_sum))  

if __name__ =='__main__':
    run__pool()