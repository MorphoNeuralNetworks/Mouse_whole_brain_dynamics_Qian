import math
import copy,os,shutil
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import pdist,cdist,squareform

def dist(A,B):
    """
    用来计算两个样本点之间的距离
    :param a: 样本点
    :param b: 样本点
    :return: 两个样本点之间的距离
    """
    temp=0
    for i in range(0,len(A)):
        temp+=(A[i]-B[i])*(A[i]-B[i])
    return math.sqrt(temp)


def CalculateDistMatrix(dataset):
    """
    计算距离矩阵
    :param dataset: 数据集
    :return: 距离矩阵
    """
    temp=pdist(dataset,metric='euclidean')
    DistMatrix = squareform(temp)

    return DistMatrix


def returnEpsCandidate(dataSet):
    """
    计算Eps候选列表
    :param dataSet: 数据集
    :return: eps候选集合
    """
    DistMatrix = CalculateDistMatrix(dataSet)
    tmp_matrix = copy.deepcopy(DistMatrix)
    for i in range(len(tmp_matrix)):
        tmp_matrix[i].sort()
    EpsCandidate = []
    for k in range(1,len(dataSet)):
        Dk = tmp_matrix[:,k]
        DkAverage = np.mean(Dk)
        EpsCandidate.append(DkAverage)
    return EpsCandidate


def returnMinptsCandidate(dataSet,DistMatrix,EpsCandidate):
    """
    计算Minpts候选列表
    :param DistMatrix: 距离矩阵
    :param EpsCandidate: Eps候选列表
    :return: Minpts候选列表
    """
    MinptsCandidate = []
    for k in range(len(EpsCandidate)):
        tmp_eps = EpsCandidate[k]
        tmp_count = np.count_nonzero(DistMatrix<=tmp_eps)
        MinptsCandidate.append(tmp_count/len(dataSet))
    return MinptsCandidate

def MutiGetDBSCAN(par):
    eps,min_samples,num,np_dataset=par[0],par[1],par[2],par[3]
    clustering = DBSCAN(eps=eps,min_samples=min_samples).fit(np_dataset)
    num_clustering = max(clustering.labels_)
    np.save('./temp/'+str(num)+'.npy',num_clustering)

def returnClusterNumberList(dataset,EpsCandidate,MinptsCandidate):
    """
    计算聚类后的类别数目 
    :param dataset: 数据集
    :param EpsCandidate: Eps候选列表
    :param MinptsCandidate: Minpts候选列表
    :return: 聚类数量列表
    """

    if os.path.exists('./temp'):
        shutil.rmtree('./temp') 
    os.mkdir('./temp')
    np_dataset = np.array(dataset)  #将dataset转换成numpy_array的形式

    ClusterNumberList = []
    par_list=[]
    for i in range(len(EpsCandidate)):
        par_list.append([EpsCandidate[i],MinptsCandidate[i],i,dataset])
    # from multiprocessing import Pool
    # cpu_worker_num = 32
    # with Pool(cpu_worker_num) as p:
    #     p.map(MutiGetDBSCAN, par_list)
    for x in par_list:
        MutiGetDBSCAN(x)
    for i in range(len(EpsCandidate)):
        tt=np.load('./temp/'+str(i)+'.npy',allow_pickle=True)
        ClusterNumberList.append(tt)
    return ClusterNumberList

if __name__ == '__main__':
    
    dataSet=np.load('temp.npy',allow_pickle=True)

    import time
    time_start = time.time()  # 记录开始时间
    EpsCandidate = returnEpsCandidate(dataSet)
    
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('EpsCandidate: '+str(time_sum))  
    
    DistMatrix = CalculateDistMatrix(dataSet)
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('DistMatrix: '+str(time_sum))  
    
    MinptsCandidate = returnMinptsCandidate(dataSet,DistMatrix,EpsCandidate)
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('MinptsCandidate: '+str(time_sum))  
    
    ClusterNumberList = returnClusterNumberList(dataSet,EpsCandidate,MinptsCandidate)
    time_end = time.time()  # 记录结束时间
    time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
    print('ClusterNumberList: '+str(time_sum))  
    '''
    '''
