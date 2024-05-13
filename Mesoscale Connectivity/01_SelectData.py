import numpy as np
import os,nrrd,csv,shutil

data_path=r'../Data/annotation_25.nrrd'
global CCFv3_model
CCFv3_model,options=nrrd.read(data_path)  # 读入 nrrd 文件
### 构建中尺度的连接矩阵
global z_half
z_half=456*25/2
global soma_info
soma_info=np.load("../Data/Soma_info.npy",allow_pickle=True).item()
soma_info={key:[soma_info[key][0],soma_info[key][1]] for key in soma_info.keys()}  
cell_type = {key: soma_info[key][0] for key in soma_info.keys()}
for key in cell_type.keys():
    if "fiber tracts" in cell_type[key]:
        cell_type[key] = "fiber_tracts"
cell_count = dict()
for key in cell_type.keys():
    if cell_type[key] not in cell_count.keys():
        cell_count[cell_type[key]] = 1
    else:
        cell_count[cell_type[key]] += 1

## 并行统计每个neuron中bouton点在不同脑区的信息
def filecalculation(par):
    root,name,wirtepath=par[0],par[1],par[2]
    empty=[]
    with open(os.path.join(root,name)) as file_object:
        contents = file_object.readlines()
        # print(len(contents))
        file_object.close()
    t=name.split('.')[0]
    while contents[0][0]=='#':# 删除注释
        del contents[0]  
    for lineid in range(0,len(contents)):
        x=contents[lineid]
        x=x.strip("\n")
        t1=x.split( )
        t1=list(map(float,t1))
        if t1[0]==0:
            continue
        if t1[0]==t1[-1]: #修正版本的数据会出现自己导向自己的情况，得处理这个问题
            print('Circle Warning!!! '+ t[0])
            continue
        if t1[1]==5:
            side=""
            if soma_info[t][1][2]>z_half:
                t1[4]=round(z_half*2-t1[4],2)
            if t1[4]<=z_half:
                side="R"
            else:
                side="L"
            tt=[round(t1[2]/25),round(t1[3]/25),round(t1[4]/25)]
            if tt[0]<528 and tt[1]<320 and tt[2]<456 and tt[2]>=0:
                empty.append([CCFv3_model[tt[0],tt[1],tt[2]],t1[2],t1[3],t1[4],side])
    # 转换为csv
    with open(wirtepath+t+".csv","w+",newline='') as f:
        csv_writer = csv.writer(f)
        for rows in empty:
            csv_writer.writerow(rows)
        f.close()

def run_pool():  # main process
    
    from multiprocessing import Pool
    cpu_worker_num = 36
    
    import time
    time_start = time.time()  # 记录开始时间

    ## 这里是含有bouton的swc数据所在的路径
    path='./bouton_swc'
    writepath='./bouton_info/'
    # 清空文件夹
    if os.path.exists(writepath):
        shutil.rmtree(writepath)  
    os.mkdir(writepath)
    # 构建参数组合
    par_list=[]
    for root,dirs,files in os.walk(path,topdown=True):
        for file in files:
            par_list.append((path,file,writepath))    
    # print(par_list[0])
    with Pool(cpu_worker_num) as p:
        p.map(filecalculation, par_list)
        
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('Select Data time: '+str(time_sum))  


if __name__ == "__main__":
    run_pool()
    # 合并区域
    path='./bouton_info/'
    bouton_region_dict=dict()
    bouton_location_dict=dict()
    for root,dirs,files in os.walk(path,topdown=True):
        for file in files:
            with open(os.path.join(path,file), newline='') as csvfile:
                contents = list(csv.reader(csvfile, delimiter=','))
            temp=[(x[0],x[4]) for x in contents]
            t=dict()
            for x in set(temp):
                t.update({x:temp.count(x)})
            tt={k:[] for k in t.keys()}
            for x in contents:
                x[1:4]=list(map(float,x[1:4]))
                tt[(x[0],x[4])].append(x[1:4])
            bouton_region_dict[file.split('.')[0]]=t
            bouton_location_dict[file.split('.')[0]]=tt
    name=""
    np.save("./"+name+"bouton_region_dict.npy",bouton_region_dict)    
    np.save("./"+name+"bouton_location_dict.npy",bouton_location_dict)  
