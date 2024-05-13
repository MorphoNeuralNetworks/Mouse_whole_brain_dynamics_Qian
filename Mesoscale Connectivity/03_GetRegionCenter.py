import numpy as np
import csv,os,random,nrrd

## 设定映射的区域范围
with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    cell_type=[str(x[0]) for x in temp]

## 寻找上层的种类
with open(r'../Data/Dataset_Structure.csv', 'r', newline='') as csvfile:
    t = csv.reader(csvfile)
    structure=np.array(list(t))
    csvfile.close()
structure_dict={}
structure_list=structure[:,3]
for x in structure:
    t=x[1].split('/')
    structure_dict[str(x[3])]=list(map(int,t[1:-1]))
structure_dict['fiber tracts']=[997, 1009]
id2name={structure_dict[x][-1]:x for x in structure_dict.keys()}
name2id={x:structure_dict[x][-1] for x in structure_dict.keys()}

## 构建cell_type与id的信息
upper_name={structure_dict[x][-1]:x for x in cell_type}
upper_id=[x for x in upper_name.keys()]
upper_dict={0:0}
## 构建向上的映射
for x in structure_dict.keys():
    path=structure_dict[x]
    mark=1
    for k in upper_id:
        if k in path:
            mark=0
            upper_dict[structure_dict[x][-1]]=k
    if mark==1:
        upper_dict[structure_dict[x][-1]]=structure_dict[x][-1]
        
## 获得脑区的中心坐标
data_path='../Data/annotation_25.nrrd'
global CCFv3_model
CCFv3_model,options=nrrd.read(data_path)  # 读入 nrrd 文件
CCFv3_model=CCFv3_model.astype(np.float64)
region_id={x:[] for x in upper_id}
x,y,z=np.shape(CCFv3_model)
for i in range(int(x)):
    if i%10==0:
        print(i)
    for j in range(int(y)):
        for k in range(int(z/2)):
            t=upper_dict[CCFv3_model[i,j,k]]
            if t in upper_id:
                region_id[t].append([i,j,k])

region_center={id2name[x]:np.mean(region_id[x],0) for x in region_id.keys()}
centres=[region_center[x].tolist() for x in cell_type]
tt=[[x[0],x[1],456-x[2]]for x in centres]
centres=centres+tt
# the method rotate the Allen 3D (x1,y1,z1) reference in the TVB 3D reference (x2,y2,z2).
# the relation between the different reference system is: x1=z2, y1=x2, z1=y2 ; x2=y1 y2=z1 z2=x1
new_centres=[[x[1],x[2],x[0]]for x in centres]
np.save("region_center.npy",new_centres)

# the method returns the tract lengths between the brain areas in the selected parcellation
def construct_tract_lengths(centres):
    len_right = len(centres) // 2
    tracts = np.zeros((len_right, len(centres)), dtype=float)
    for inj in range(len_right):
        center_inj = centres[inj]
        for targ in range(len_right):
            targ_r = centres[targ]
            targ_l = centres[targ + len_right]
            tracts[inj, targ] = np.sqrt(
                (center_inj[0] - targ_r[0]) ** 2 + (center_inj[1] - targ_r[1]) ** 2 + (center_inj[2] - targ_r[2]) ** 2)
            tracts[inj, targ + len_right] = np.sqrt(
                (center_inj[0] - targ_l[0]) ** 2 + (center_inj[1] - targ_l[1]) ** 2 + (center_inj[2] - targ_l[2]) ** 2)
    # Save the complete matrix (both left and right inj):
    first_quarter = tracts[:, :(tracts.shape[1] // 2)]
    second_quarter = tracts[:, (tracts.shape[1] // 2):]
    tracts_down = np.concatenate((second_quarter, first_quarter), axis=1)
    tracts = np.concatenate((tracts, tracts_down), axis=0)
    return tracts.T  # transpose because TVB convention requires SC[target, source]!
tract_lengths=construct_tract_lengths(centres)
np.save("tract_lengths.npy",tract_lengths)
