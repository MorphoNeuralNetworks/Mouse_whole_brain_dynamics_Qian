import numpy as np
import csv,os,shutil

global actual_sw,actual_sp
actual_sw,actual_sp=60,2 # unit second

def FCDGet(par):
    (coupling_v,noise_v)=par
    filepath='./BOLD_single/'

    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    bold_time=np.load(filepath+"bold_time"+file_end+".npy",allow_pickle=True)
    bold_data=np.load(filepath+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data=bold_data[:,0,:,0].T
    
    # start time and end time
    padding_s,padding_e=60,660
    
    bold_time=np.array([(x+1)*1000 for x in range(0,len(bold_time)-padding_s*2)])
    bold_data=bold_data[:,padding_s:padding_e]
    
    node_len,time_len=bold_data.shape
    
    steps=int(np.floor((time_len-actual_sw)/actual_sp))
    
    fc=np.zeros((steps,node_len,node_len))
    for i in range(steps):
        st_p,end_p=actual_sp*i,actual_sp*i+actual_sw
        cut_data=bold_data[:,st_p:end_p]
        cor_t=np.corrcoef(cut_data)
        fc[i,:,:]=cor_t
    temp=[]
    for i in range(steps):
        fc_t=fc[i,:]
        triangular = np.triu_indices(len(fc_t), 1)
        tt=fc_t[triangular]
        temp.append(tt)
    temp=np.array(temp)
    fcd=np.corrcoef(temp)
    wirtepath="./FCD_single/"
    np.save(wirtepath+"FCD"+file_end+".npy",fcd)
    
def run__pool():  # main process
    from multiprocessing import Pool
    cpu_worker_num = 36 #设定并行的数量
    par_list=[]
    par_list=[]
    
    # set the parameters range
    for i in range(0,100):
        for j in range(0,100):
            par_list.append((i*0.001+0.2,j*0.000001+0.000001))
    
    writepath='./FCD_single/'
    if os.path.exists(writepath):
        shutil.rmtree(writepath)  
    os.mkdir(writepath) 
        
    import time
    time_start = time.time()  # 记录开始时间
        
    with Pool(cpu_worker_num) as p:
        p.map(FCDGet, par_list)
    # for x in par_list:
    #     FCDGet(x)
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('Get FCD time: '+str(time_sum))  

if __name__ =='__main__':
    run__pool()

