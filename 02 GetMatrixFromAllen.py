#### Use allen's mesoscale experimental data to obtain connections between specified regions
import numpy as np
import csv,json

## CCFv3 structure
with open('./Data/tree.json','r',encoding='utf_8_sig')as fp:
    json_data = json.load(fp)
    fp.close()
Celltype2Id={x['acronym']:x['id'] for x in json_data}
Fullname2Type={x['name']:x['acronym'] for x in json_data}
Id2Celltype={Celltype2Id[key]:key for key in Celltype2Id.keys()}

## Set regions interested
with open('./Data/TVMB_Regions.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Region=np.array(temp)
TVMB_Re=[]
for x in TVMB_Region:
    if x[0].replace("_"," ") not in Fullname2Type.keys():
        print(x[0])
    else:
        TVMB_Re.append(Fullname2Type[x[0].replace("_"," ")])

## From a divided region to a larger region
# select_region=[]
# Celltype2Select=dict()
# for x in json_data:
#     mark=0
#     for k in select_region:
#         if Celltype2Id[k] in x['structure_id_path']:
#             Celltype2Select[x['acronym']]=k
#             mark=1
#             break
#     if mark==0:
#         Celltype2Select[x['acronym']]=x['acronym']

def Delete3std(x): ## delete value out of 2*std
    mean=np.mean(x)
    std=np.std(x)
    temp=[]
    for t in x:
        if t>=mean-2*std and t<=mean+2*std:
            temp.append(t)
    return np.mean(temp)

## Get all Allen experimental data
with open('./Data/mouse_projection_data_sets.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=','))
    del temp[0]
    dataset=np.array(temp)

all_exp=np.zeros((len(TVMB_Re)*2,len(TVMB_Re)*2))
injection_density=dict()
for i in range(0,len(TVMB_Re)):
    print(i)
    check=[]
    i_r=TVMB_Re[i] # set regions
    i_r_id=Celltype2Id[i_r]

    exp_list=[]
    if i_r=="CENT" or i_r=="CUL" or i_r=="AN": # some special regions without experiments
        for x in dataset:
            if i_r in x[2]: # Contain subregions
                exp_list.append(str(x[0])) 
    else:
        for x in dataset:
            # if x[2]==i_r: # fouce on primary injection structure
            #     exp_list.append(str(x[0]))  
            tt=x[4][1:-1].split(',')
            if i_r in tt: # fouce on secondary injection structure
                exp_list.append(str(x[0]))  

    i_r=TVMB_Re[i] # set region 
    exp_r=[]
    for exp in exp_list:
        with open('./Data/meso_projection/experiment_'+exp+'.csv', newline='') as csvfile:
            temp = list(csv.reader(csvfile, delimiter=','))
        del temp[0]
        # Right hemisphere; Injection voxels > 50; Projection volume > 2
        mark=0
        for tt in temp:
            if int(tt[2])==i_r_id and tt[4]=='t' and tt[3]=='2':
                check.append([tt[6],tt[13]])
                if float(tt[6])<200 or float(tt[13])<0.03125:
                # if float(tt[6])<0 or float(tt[13])<0:
                    # print(exp)
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

    # find all suitable experiments
    for j in range(0,len(TVMB_Re)):
        p_r=TVMB_Re[j]
        p_r_id=Celltype2Id[p_r]
        exp_t_l=[]
        exp_t_r=[]
        # get connectivity
        for x in exp_r:
            if int(x[2])==p_r_id:
                if x[3]=='1':
                    exp_t_l.append(float(x[6])/float(x[5])/injection_density[x[1]])
                if x[3]=='2':
                    exp_t_r.append(float(x[6])/float(x[5])/injection_density[x[1]])
       
        if len(exp_t_r)!=0:
            all_exp[i,j]=Delete3std(exp_t_r)
        if len(exp_t_l)!=0:
            all_exp[i,j+len(TVMB_Re)]=Delete3std(exp_t_l)

# mirror
all_exp[len(TVMB_Re):,len(TVMB_Re):]=all_exp[0:len(TVMB_Re),0:len(TVMB_Re)]
all_exp[len(TVMB_Re):,0:len(TVMB_Re)]=all_exp[0:len(TVMB_Re),len(TVMB_Re):]

np.save("TVMB_connectivity_Allen.npy",all_exp)    