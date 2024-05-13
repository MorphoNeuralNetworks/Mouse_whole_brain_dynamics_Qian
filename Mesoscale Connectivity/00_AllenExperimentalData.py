import numpy as np
import csv,json,nrrd,os

## 统计每个区域内的大小
data_path=r'../Data/annotation_25.nrrd'
CCFv3_model,options=nrrd.read(data_path)  # 读入 nrrd 文件
pixel_num=[0]*(np.max(CCFv3_model)+1)
for i in range(0,np.shape(CCFv3_model)[0]):
    if i%10==0:
        print(i)
    for j in range(0,np.shape(CCFv3_model)[1]):
        # if j%100==0:
        #     print(j)
        for k in range(0,np.shape(CCFv3_model)[2]):
            if CCFv3_model[i,j,k]!=0:
                pixel_num[CCFv3_model[i,j,k]]+=1
pixel_dict=dict()
for i in range(0,len(pixel_num)):
    if pixel_num[i]!=0:
        pixel_dict[i]=pixel_num[i]
np.save("pixel_num_dict.npy",pixel_dict)

## 构建CCFv3结构
with open('../Data/tree.json','r',encoding='utf_8_sig')as fp:
    json_data = json.load(fp)
    fp.close()
Celltype2Id={x['acronym']:x['id'] for x in json_data}
Fullname2Type={x['name']:x['acronym'] for x in json_data}
Id2Celltype={Celltype2Id[key]:key for key in Celltype2Id.keys()}

## 设定关注的脑区 neuron>5/10/20的三种情况
with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]

def Delete3std(x): #删除2std以外的数据
    mean=np.mean(x)
    std=np.std(x)
    temp=[]
    for t in x:
        if t>=mean-2*std and t<=mean+2*std:
            temp.append(t)
    return np.mean(temp)

with open('../Data/mouse_projection_data_sets.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=','))
    del temp[0]
    dataset=np.array(temp)

all_exp=np.zeros((len(TVMB_Re)*2,len(TVMB_Re)*2))
injection_density=dict()
for i in range(0,len(TVMB_Re)):
    print(i)
    check=[]
    i_r=TVMB_Re[i] #设置筛选区域
    i_r_id=Celltype2Id[i_r]

    exp_list=[]
    mark=0
    for x in dataset:
        # if x[2]==i_r: # 只看一级结构匹配度
        #     exp_list.append(str(x[0]))  
        tt=x[4][1:-1].split(',')
        if i_r in tt: # 看二级结构匹配度
            exp_list.append(str(x[0]))  
            mark=1
    if mark==0:
        print(TVMB_Re[i])

    i_r=TVMB_Re[i] #设置筛选区域
    exp_r=[]
    for exp in exp_list:
        with open('../Data/meso_projection/experiment_'+exp+'.csv', newline='') as csvfile: #这是Allen的2062个注射实验
            temp = list(csv.reader(csvfile, delimiter=','))
        del temp[0]
        # 只看右半球 只看注射面积大于50
        mark=0
        for tt in temp:
            if int(tt[2])==i_r_id and tt[4]=='t' and tt[3]=='2':
                check.append([tt[6],tt[13]])
                if float(tt[6])<200 or float(tt[13])<0.03125:
                    break
                else:
                    injection_density[tt[1]]=float(tt[6])/float(tt[5])
                    mark=1

        if mark==1:
            tempi=[]
            for tt in temp:
                if tt[4]=='f' and tt[3]!='3':
                    tempi.append(tt)    
            exp_r=exp_r+tempi
 
    exp_r=np.array(exp_r)

    # 找到目标区域符合的实验
    for j in range(0,len(TVMB_Re)):
        p_r=TVMB_Re[j] #设置筛选区域
        p_r_id=Celltype2Id[p_r]
        exp_t_l=[]
        exp_t_r=[]
        for x in exp_r:
            if int(x[2])==p_r_id:
                if x[3]=='1':
                    exp_t_l.append(float(x[6])/float(x[5])/injection_density[x[1]])
                if x[3]=='2':
                    # exp_t.append(float(x[14]))
                    exp_t_r.append(float(x[6])/float(x[5])/injection_density[x[1]])
       
        if len(exp_t_r)!=0:
            all_exp[i,j]=Delete3std(exp_t_r)
        if len(exp_t_l)!=0:
            all_exp[i,j+len(TVMB_Re)]=Delete3std(exp_t_l)

# transpose because TVB convention requires SC[target, source]!
all_exp[0:len(TVMB_Re),0:len(TVMB_Re)]=all_exp[0:len(TVMB_Re),0:len(TVMB_Re)].T
all_exp[0:len(TVMB_Re),len(TVMB_Re):]=all_exp[0:len(TVMB_Re),len(TVMB_Re):].T
# 镜像对称
all_exp[len(TVMB_Re):,len(TVMB_Re):]=all_exp[0:len(TVMB_Re),0:len(TVMB_Re)]
all_exp[len(TVMB_Re):,0:len(TVMB_Re)]=all_exp[0:len(TVMB_Re),len(TVMB_Re):]

np.save("./Allen_connectivity.npy",all_exp)    
