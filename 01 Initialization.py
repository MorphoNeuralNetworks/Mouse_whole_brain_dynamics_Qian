#### Shared data to be used
import numpy as np
import nrrd

## get size of different brain regions
data_path=r'./Data/annotation_25.nrrd'
CCFv3_model,options=nrrd.read(data_path)
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
