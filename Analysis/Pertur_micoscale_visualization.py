import numpy as np
import csv

with open('./Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]
area_color=np.load(r'D:\QPH\Data\Other_Infomation\color_network.npy', allow_pickle=True).item()

# names=["","scale_0.5_all_","scale_0.5_axon_","scale_0.5_dendrite_","prune_0.5_all_","prune_0.5_axon_","prune_0.5_dendrite_","delete_0.5_all_"]
names=["","scale_0.5_all_","prune_0.5_all_","delete_0.5_all_"]
for name in names[0:1]:
    # name=""
    connection=np.load("./Connectivity/"+name+"SingleCell_connectivity_10.npy",allow_pickle=True)
    # connection=connection[0:len(TVMB_Re),:]
    # for i in range(0,len(TVMB_Re)): connection[i,i]=0 
    connection=np.log(connection+1)
    # connection=(connection-np.min(connection))/(np.max(connection)-np.min(connection))
    '''
    ## heatmap
    from matplotlib import pyplot as plt
    import seaborn as sns
    plt.close("all")
    color_list=[area_color[x] for x in TVMB_Re]+[area_color[x] for x in TVMB_Re]
    plot=sns.clustermap(connection[0:len(TVMB_Re),:],figsize=(10,6),xticklabels=TVMB_Re,yticklabels=TVMB_Re,
                        row_cluster=False,col_cluster=False,annot=False,
                        row_colors=color_list,col_colors=color_list,square=True,
                        cbar_kws={'shrink':0.3,'orientation': 'horizontal'},
                        annot_kws={'fontsize':10})
    plt.tight_layout()
    plt.savefig("./Connectivity/"+name+'heatmap.jpg', dpi=300)
    '''
    
    ##igraph
    import igraph as ig
    import cairo
    import cv2
    import math
    soma_info=np.load(r'D:\QPH\Data_new\Other_Infomation\Soma_info.npy',allow_pickle=True).item()  
    cell_type=[soma_info[key][0] for key in soma_info.keys()]
    cell_count={}
    for item in cell_type:
        cell_count.update({item:cell_type.count(item)})
    
    region_center=np.load("region_center.npy",allow_pickle=True)
    
    node_list=[x+"_R" for x in TVMB_Re]+[x+"_L" for x in TVMB_Re]
    node_color=[area_color[x] for x in TVMB_Re]+[area_color[x] for x in TVMB_Re]
    node_size=[np.log(cell_count[x])*6 for x in TVMB_Re]+[np.log(cell_count[x])*6 for x in TVMB_Re]
    
    edge_list=[]
    edge_weight=[]
    edge_width=[]
    edge_color=[]
    for i in range(0,len(TVMB_Re)):
        for j in range(0,len(node_list)):
            if connection[i,j]!=0:
                edge_list.append((node_list[i],node_list[j]))
                edge_weight.append(connection[i,j])
                edge_width.append(connection[i,j])
                edge_color.append(area_color[TVMB_Re[i]])
    
    node_list.append('empty point_up')
    node_list.append('empty point_down')
    node_local=region_center[:,[1,2]].tolist()
    node_local.append((0,0))
    # node_local.append((528,456))
    node_local.append((456,528))
    node_size.append(0)
    node_size.append(0)
    
    g = ig.Graph(directed=True)
    
    g.add_vertices(node_list)
    g.add_edges(edge_list)
    # g.es['weight']=edge_weight
    g.es['width']=edge_width
    g.es['color']=edge_color
    g.es['arrow_size']=0.3
    
    g.vs['color'] = node_color
    g.vs['label'] = node_list
    g.vs['location']=node_local
    g.vs['size']=node_size
    g.vs['label_size']=[15]*len(TVMB_Re)*2+[0,0]
    print(ig.summary(g))
    
    # 528 320 456
    # out=ig.plot(g,layout="circle",bbox = (1500,1500))
    # out.save("./Connectivity/"+name+'NetworkCircle.png')
    
    # 528 320 456
    # out=ig.plot(g,layout=g.vs['location'],bbox = (528*2,456*2))
    out=ig.plot(g,layout=g.vs['location'],bbox = (456*2,528*2))
    out.save("./Connectivity/"+name+'Network.png')
    
    img1=cv2.imread('background.png')
    img2=cv2.imread("./Connectivity/"+name+'Network.png')
    res = img1
    for i in range(0,np.shape(img2)[0]):
        for j in range(0,np.shape(img2)[1]):
            if np.sum(img2[i,j])!=255*3:
                res[i,j]=img2[i,j]
    # 保存
    # cv2.imencode('.png', res)[1].tofile("./Connectivity/"+name+"BrainNetwork.png")
