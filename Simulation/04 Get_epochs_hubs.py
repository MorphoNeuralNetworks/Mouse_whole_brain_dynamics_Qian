import numpy as np
import csv,os,shutil
import KANN_Dbscan as KD
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.spatial.distance import pdist
from sklearn.cluster import DBSCAN
from sklearn.manifold import SpectralEmbedding
import igraph as ig

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def GetDbscanPar(data):
    EpsCandidate = KD.returnEpsCandidate(data)
    DistMatrix = KD.CalculateDistMatrix(data)
    MinptsCandidate = KD.returnMinptsCandidate(data,DistMatrix,EpsCandidate)
    ClusterNumberList = KD.returnClusterNumberList(data,EpsCandidate,MinptsCandidate)
    for i in range(0,len(ClusterNumberList)-3):
        if ClusterNumberList[i]==ClusterNumberList[i+1] and ClusterNumberList[i+1]==ClusterNumberList[i+2]:
            if ClusterNumberList[i]!=0:
                count=0
                while(ClusterNumberList[i+count]==ClusterNumberList[i+1+count]):
                    count+=1
                best_eps=EpsCandidate[i+count]
                best_min_samples=MinptsCandidate[i+count]
                break
            else:
                print("Warning cluster==0")
                best_eps=EpsCandidate[i+2]
                best_min_samples=MinptsCandidate[i+2]
                break
    return best_eps,best_min_samples

def spectral_dbscan(fcd, n_dim=2):
    fcd = fcd - fcd.min()
    se = SpectralEmbedding(n_dim, affinity="precomputed")
    xi = se.fit_transform(fcd)
    pd = pdist(xi)
    best_eps,best_min_samples=GetDbscanPar(xi)
    db = DBSCAN(eps=best_eps,min_samples=best_min_samples).fit(xi)
    return xi.T, db.labels_

def compute_radii(xi, centered=False):
    if centered:
        xi = xi.copy() - xi.mean(axis=1).reshape((len(xi), 1))
    radii = np.sqrt(np.sum(xi ** 2, axis=0))
    return radii

def spectral_embedding(fcd):
    xi, _ = spectral_dbscan(fcd, 2)
    xir = compute_radii(xi, True)
    xir_sorted = np.sort(xir)
    xir_cutoff = 0.5 * xir_sorted[-1]
    return xir, xir_cutoff

def epochs_interval(xir, xir_cutoff, sp, sw):
        # Calculate the starting point and the ending point of each epoch of stability
        # sp=spanning, sw=sliding window
        epochs_dict = {}  # here the starting and the ending point will be stored
        thresholds = np.where(xir < xir_cutoff)
        tt = 0
        ep = 0
        while (tt + 2) < len(thresholds[0]):
            epochs_dict[ep] = [thresholds[0][tt]]  # starting point of epoch ep
            while ((tt + 2) != len(thresholds[0])) & (thresholds[0][tt + 1] == thresholds[0][tt] + 1):
                # until the vector is not finish and until each entries +1 is equal to the next one
                tt += 1
            epochs_dict[ep].append(thresholds[0][tt])
            tt += 1
            ep += 1

        epochs_extremes = np.zeros((len(epochs_dict), 2), dtype=float)
        for ep in range(len(epochs_dict)):
            epochs_extremes[ep, 0] = epochs_dict[ep][0] * sp
            epochs_extremes[ep, 1] = epochs_dict[ep][1] * sp + sw
        return epochs_extremes

def nodal_eff(g):
    weights = g.es["weight"][:]
    sp = (1.0 / np.array(g.shortest_paths_dijkstra(weights=weights)))
    np.fill_diagonal(sp,0)
    N=sp.shape[0]
    ne= (1.0/(N-1)) * np.apply_along_axis(sum,0,sp)

    return np.mean(ne)

def shortest_path(g):
    cc=np.array(g.shortest_paths_dijkstra(weights=g.es["weight"][:]))
    return np.sum(cc)/(len(cc)*(len(cc)-1))

global actual_sw,actual_sp,num_eig
actual_sw,actual_sp=60,2 # unit second
num_eig = 2  # number of the eigenvector that will be extracted
single_par,allen_par=np.load("best_par.npy",allow_pickle=True)
for i in range(len(allen_par)):
    if allen_par[i,0]<0.084:
        allen_par[i,0]=0.084

with open(r'..\Data\SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0]) for x in temp]
full_TVMB_Re=[x+"_R" for x in TVMB_Re]+[x+"_L" for x in TVMB_Re]
full_TVMB_Re=np.array(full_TVMB_Re)

result={}
# 对每个实验计算稳态的epoch
for exp_id in range(0,20):
    par=single_par[exp_id]
    writepath="./kdeplot_single/"
    filepath='./FCD_single/'
    coupling_v,noise_v=par[0],par[1]
    file_end="_w_"+str(round(coupling_v,4))+"_n_"+str(round(noise_v,6))
    fcd=np.load(filepath+"FCD"+file_end+".npy",allow_pickle=True)
    bold_data=np.load('./BOLD_single_detail/'+"bold_data"+file_end+".npy",allow_pickle=True)
    bold_data=bold_data[:,0,:,0].T
    
    eigvect_dict = []  # holds eigenvectors of the fcs calculated over the epochs, key1=mode, key2=var, key3=numb ep
    eigval_dict = []  # holds eigenvalues of the fcs calculated over the epochs, key1=mode, key2=var, key3=numb ep
    fcd_matrix=fcd
    
    [xir, xir_cutoff] = spectral_embedding(fcd_matrix)
    plt.plot(xir)
    plt.plot([xir_cutoff]*len(xir))
    
    epochs_extremes = epochs_interval(xir, xir_cutoff, actual_sp, actual_sw)
    for i in range(len(epochs_extremes)-1,-1,-1):
        if epochs_extremes[i][1]-epochs_extremes[i][0]<(10+actual_sw):
            epochs_extremes=np.delete(epochs_extremes,i,axis=0)
    xir_mark=xir
    
    epoch_len=[x[1]-x[0] for x in epochs_extremes]
    tt=epoch_len.index(np.max(epoch_len))
    st_p,end_p=int(epochs_extremes[tt][0]), int(epochs_extremes[tt][1]) + 1
    
    cut_data=bold_data[:,st_p:end_p]
    fc = np.corrcoef(cut_data)  # calculate fc over the epoch of stability
    
    eigval_matrix, eigvect_matrix = linalg.eig(fc)
    eigval_matrix = np.real(eigval_matrix)
    eigvect_matrix = np.real(eigvect_matrix)
    eigval_matrix = eigval_matrix / np.sum(np.abs(eigval_matrix))  # normalize eigenvalues to [0 and 1)
    for en in range(num_eig):
        index = np.argmax(eigval_matrix)
        eigvect_dict.append(abs(eigvect_matrix[:, index]))
        eigval_dict.append(eigval_matrix[index])
        eigval_matrix[index] = 0
    for tt in range(0,num_eig):
        eig=eigvect_dict[tt]
        print("Eigvector "+str(tt)+", Weight "+str(eigval_dict[tt]))
        cc=full_TVMB_Re[np.argsort(-eig)]
        dd=eig[np.argsort(-eig)]
        result[(exp_id,tt)]=np.array([[cc[i],dd[i]] for i in range(len(cc))])

    conn=[]
    for i in range(0,len(fc)):
        for j in range(i+1,len(fc)):
            if fc[i,j]!=0:
                conn.append((full_TVMB_Re[i],full_TVMB_Re[j],fc[i,j]))
    data=np.array(conn)
    
    node_list=list(set(data[:,0].tolist()+data[:,1].tolist()))
    edge_list=[(str(x[0]),str(x[1])) for x in data]
    edge_list=list(set(edge_list))
    weight_dict=dict()
    for x in data:
        if (str(x[0]),str(x[1])) not in weight_dict.keys():
            weight_dict[(str(x[0]),str(x[1]))]=1/np.absolute(float(x[2]))
        else:
            print("error")
            weight_dict[(str(x[0]),str(x[1]))]+=float(x[2])
    edge_weight=[weight_dict[x] for x in edge_list]
       
    g = ig.Graph(directed=False)
    g.add_vertices(node_list)
    g.add_edges(edge_list)
    g.es['weight']=edge_weight
    result[exp_id]=[g.transitivity_avglocal_undirected(weights=g.es["weight"][:]),
                           shortest_path(g),
                           nodal_eff(g)]

## 可视化划分的epoch
# 取反向区间
epochs_temp=np.array([[x[0]/actual_sp,(x[1]-actual_sw)/actual_sp]for x in epochs_extremes])
epochs_reverse=[]
if epochs_temp[0,0]!=0:
    epochs_reverse.append([0,int(epochs_temp[0,0])])
for i in range(1,len(epochs_temp)):
    epochs_reverse.append([int(epochs_temp[i-1,1]),int(epochs_temp[i,0])])
if epochs_extremes[-1,1]!=len(fcd):
    epochs_reverse.append([int(epochs_extremes[-1,1]),len(fcd)])
epochs_reverse=np.array(epochs_reverse)
fcd_segmented = fcd.copy()

for i in range(0,len(epochs_reverse)):
    st,ed=int(epochs_reverse[i,0]),int(epochs_reverse[i,1])
    fcd_segmented[st:ed,:] = 1.1
    fcd_segmented[:,st:ed] = 1.1

plt.close("all")
plt.figure(figsize=(9.5,4))
plt.subplot(1,2,1)
cs=plt.imshow(fcd,cmap='jet',vmin=0,vmax=1)
axcb=plt.colorbar(ticks=[0, 0.5, 1])
cs.set_clim(0, 1.0)
plt.xlabel(r'Time $t_j$ (s)', fontsize=12)
plt.ylabel(r'Time $t_i$ (s)', fontsize=12)
plt.title('FCD', fontsize=14)

plt.subplot(1,2,2)
cs=plt.imshow(fcd_segmented,cmap='jet',vmin=0)
axcb=plt.colorbar(ticks=[0, 0.5, 1])
cs.set_clim(0, 1.0)
plt.xlabel(r'Time $t_j$ (s)', fontsize=12)
plt.ylabel(r'Time $t_i$ (s)', fontsize=12)
plt.title('FCD segmented', fontsize=14)
plt.tight_layout()

## 计算epoch分段的eigvalue和eigvector
for ep in range(0, epochs_extremes.shape[0]):
    if epochs_extremes[0,0]==0 and ep==0:
        continue
    eigvect_dict[ep] = []
    eigval_dict[ep] = []
    
    st_p,end_p=int(epochs_extremes[ep][0]), int(epochs_extremes[ep][1]) + 1
    cut_data=bold_data[:,st_p:end_p]
    fc = np.corrcoef(cut_data)  # calculate fc over the epoch of stability
    
    eigval_matrix, eigvect_matrix = linalg.eig(fc)
    eigval_matrix = np.real(eigval_matrix)
    eigvect_matrix = np.real(eigvect_matrix)
    eigval_matrix = eigval_matrix / np.sum(np.abs(eigval_matrix))  # normalize eigenvalues to [0 and 1)
    for en in range(num_eig):
        index = np.argmax(eigval_matrix)
        eigvect_dict[ep].append(abs(eigvect_matrix[:, index]))
        eigval_dict[ep].append(eigval_matrix[index])
        eigval_matrix[index] = 0

## hubs的可视化
x_num,y_num=len(eigval_dict),num_eig
fig, axes = plt.subplots(x_num,y_num)

j=0
import csv
with open('../Data/SEU_Regions_10.csv', newline='') as csvfile:
    temp = list(csv.reader(csvfile, delimiter=' '))
    TVMB_Re=[str(x[0])+'_R' for x in temp]+[str(x[0])+'_L' for x in temp]
    TVMB_Re_Raw=[str(x[0]) for x in temp]
    
with open(r'..\Data\Dataset_Structure.csv', 'r', newline='') as csvfile:
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
for kk in range(0,len(epochs_extremes)):
    if epochs_extremes[0,0]==0 and kk==0:
        continue
    for tt in range(0,num_eig):
        CCF_model=np.load("../Data/CCF_model.npy",allow_pickle=True)
        CCF_model=CCF_model[:,:,::-1]
        template=np.load("template.npy",allow_pickle=True)
        eig=eigvect_dict[kk][tt]
        print("Epoch "+str(kk)+", Eigvector "+str(tt)+", Weight "+str(eigval_dict[kk][tt]))
        for i in range(len(TVMB_Re)):
            if eig[i]>0.23:
                print([TVMB_Re[i],eig[i]])
        print("\n")

        Vol=np.ones(np.shape(CCF_model))*-1
        
        for i in range(len(TVMB_Re)):
            if i<len(TVMB_Re)/2:
                Vol[np.where(CCF_model[0:57,:,:]==name2id[TVMB_Re_Raw[i]])]=eig[i]
                CCF_model[np.where(CCF_model[0:57,:,:]==name2id[TVMB_Re_Raw[i]])]=eig[i]
            else:
                Vol[np.where(CCF_model==name2id[TVMB_Re_Raw[i-len(TVMB_Re_Raw)]])]=eig[i]
    
        Vol = np.ma.masked_where(Vol < 0, Vol)
        slice_idy=68
        im1 = axes[int(np.floor(j/num_eig)),j%num_eig].imshow((template[:,slice_idy,:].T)[::-1], cmap='gray', vmin=template.min(), vmax=template.max())
        cax = axes[int(np.floor(j/num_eig)),j%num_eig].imshow((Vol[:,slice_idy,:].T)[::-1], cmap='YlOrRd', alpha=1, vmin=0.15, vmax=np.amax(eig))
    
        # axes[int(np.floor(j/3)),j%3].axis('off')
        # axes[int(np.floor(j/3)),j%3].set_title(conn_measure.title)
        # divider = make_axes_locatable(axes[int(np.floor(j/3)),j%3])
        # cax1 = divider.append_axes("right", size="5%", pad=0.05)
        # axcb=plt.colorbar(cax,cax1,ticks=[0.15,np.amax(eig)],orientation='vertical')
        # axcb.set_ticklabels(['0.15',str(np.round(np.amax(eig),2))]) 
        # axcb.set_label('Eigenvector components')
        j=j+1


