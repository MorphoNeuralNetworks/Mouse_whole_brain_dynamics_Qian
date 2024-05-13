import csv,os,json
import numpy as np
import nrrd
from scipy import stats
import seaborn as sns

## 选择的区域
with open('../Data/Selected_Regions.csv', 'r', newline='',encoding='utf8') as csvfile:
    t = csv.reader(csvfile)
    t=list(t)
    select_region=[x[2] for x in t[1:]]
    csvfile.close()

## 构建CCFv3结构
with open('../Data/tree.json','r',encoding='utf_8_sig')as fp:
    json_data = json.load(fp)
    fp.close()
Celltype2Id={x['acronym']:x['id'] for x in json_data}
Id2Celltype={Celltype2Id[key]:key for key in Celltype2Id.keys()}

Celltype2Select=dict()
for x in json_data:
    mark=0
    for k in select_region:
        if Celltype2Id[k] in x['structure_id_path']:
            Celltype2Select[x['acronym']]=k
            mark=1
            break
    if mark==0:
        Celltype2Select[x['acronym']]=x['acronym']
names=[""]
for name in names:
    # 得到所有区域内bouton的数量
    bouton_region_dict=np.load("../"+name+"bouton_region_dict.npy",allow_pickle=True).item()
    bouton_region_select=dict()
    for key in bouton_region_dict.keys():
        temp=bouton_region_dict[key]
        if ('0','L') in temp.keys():
            temp.pop(('0','L'))
        if ('0','R') in temp.keys():
            temp.pop(('0','R'))
        temp={(Id2Celltype[int(key[0])],key[1]):temp[key] for key in temp.keys()}
        tt=dict()
        for k in temp.keys():
            kk=(Celltype2Select[k[0]],k[1])
            if kk not in tt.keys():
                tt[kk]=0
            tt[kk]=tt[kk]+temp[k]
        bouton_region_select[key]=tt
    
    # 得到pixel的数量
    pixel_num_dict=np.load("pixel_num_dict.npy",allow_pickle=True).item()
    pixel_num_dict={Id2Celltype[int(key)]:pixel_num_dict[key] for key in pixel_num_dict.keys()}
    pixel_num_select=dict()
    tt=dict()
    for k in pixel_num_dict.keys():
        if Celltype2Select[k] not in pixel_num_select.keys():
            pixel_num_select[Celltype2Select[k]]=0
        pixel_num_select[Celltype2Select[k]]=pixel_num_select[Celltype2Select[k]]+pixel_num_dict[k]
    
    # 得到neuron的数量
    soma_info=np.load("../Data/Soma_info.npy",allow_pickle=True).item()
    soma_info={key:soma_info[key][0] for key in soma_info.keys()}  
    tt=list(set([soma_info[key] for key in soma_info.keys()]))
    soma_num={x:0 for x in tt}
    for key in soma_info.keys():
        soma_num[soma_info[key]]=soma_num[soma_info[key]]+1
    cell_type=dict()
    for key in soma_info.keys():
        if soma_info[key] not in cell_type.keys():
            cell_type[soma_info[key]]=[]
        cell_type[soma_info[key]].append(key)
    
    # (bouton num / pixel in target region) / (soma num / pixel in source region)
    with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
        temp = list(csv.reader(csvfile, delimiter=' '))
        TVMB_Re=[str(x[0]) for x in temp]
    
    connection=np.zeros((len(TVMB_Re)*2,len(TVMB_Re)*2))
    for i in range(0,len(TVMB_Re)):
    # for i in range(22,23):
        s_r=TVMB_Re[i]
        for j in range(0,len(TVMB_Re)):
        # for j in range(22,23):
            p_r=TVMB_Re[j]
            if s_r in cell_type.keys():
                bouton_num_right,bouton_num_left=0,0
                for neuron in cell_type[TVMB_Re[i]]: 
                    if (p_r,'R') in bouton_region_select[neuron].keys(): # right
                        bouton_num_right+=bouton_region_select[neuron][(p_r,'R')]
                    if (p_r,'L') in bouton_region_select[neuron].keys(): # left                   
                        bouton_num_left+=bouton_region_select[neuron][(p_r,'L')]
            connection[i,j]=(bouton_num_right/pixel_num_select[p_r])/(soma_num[s_r]/pixel_num_select[s_r])
            connection[i,j+len(TVMB_Re)]=(bouton_num_left/pixel_num_select[p_r])/(soma_num[s_r]/pixel_num_select[s_r])
    # transpose because TVB convention requires SC[target, source]!
    connection[0:len(TVMB_Re),0:len(TVMB_Re)]=connection[0:len(TVMB_Re),0:len(TVMB_Re)].T
    connection[0:len(TVMB_Re),len(TVMB_Re):]=connection[0:len(TVMB_Re),len(TVMB_Re):].T
    # 镜像对称
    connection[len(TVMB_Re):,len(TVMB_Re):]=connection[0:len(TVMB_Re),0:len(TVMB_Re)]
    connection[len(TVMB_Re):,0:len(TVMB_Re)]=connection[0:len(TVMB_Re),len(TVMB_Re):]

    np.save("./"+name+"SingleCell_connectivity.npy",connection)  


