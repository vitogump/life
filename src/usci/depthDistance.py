# -*- coding: UTF-8 -*-
'''
Created on 2019��7��30��

@author: liurui
'''

from optparse import OptionParser
import os,re
import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline as spline
from scipy.spatial.distance import pdist

parser = OptionParser()
parser.add_option("-d", "--flexDir", dest="flexDir", help="")
parser.add_option("-t", "--distype", dest="distype",default=None, help="")
parser.add_option("-s", "--selectbinfiles", dest="selectbinfiles")

parser.add_option("-o", "--outputprefix", dest="outputpre")

(options, args) = parser.parse_args()
flexDir=options.flexDir
flexfilelist=[];flexdatalist=[]
f=open(options.selectbinfiles,'r')
selectedidx=[int(e.strip()) for e in f.readlines()];f.close()

def tryint(s):
    try:
        return int(s)
    except ValueError:
        return s
def str2int(v_str):
    return [tryint(sub_str) for sub_str in re.split('[^0-9]+',v_str)]

def cosine_distance(matrix1,matrix2):
    #this code may have some problem in python 3
    matrix1_matrix2 = np.dot(matrix1, matrix2.transpose())
    matrix1_norm = np.sqrt(np.multiply(matrix1, matrix1).sum(axis=1))
    matrix1_norm = matrix1_norm[:, np.newaxis]
    
    matrix2_norm = np.sqrt(np.multiply(matrix2, matrix2).sum(axis=1))
    matrix2_norm = matrix2_norm[:, np.newaxis]
    
    cosine_distance = np.divide(matrix1_matrix2, np.dot(matrix1_norm, matrix2_norm.transpose()))
    return cosine_distance

mds_scale=[["0" for i in range(75)] for j in range(75)]
mds_correct=[["0" for i in range(75)] for j in range(75)]

ofs=open(options.outputpre+str(options.distype)+"scaledDepth","w")
ofc=open(options.outputpre+str(options.distype)+"GCcorrectedDepth","w")
if __name__ == '__main__':
    #read file
    for elem in os.listdir(path=flexDir):
        path = flexDir + os.sep + elem
        if ( os.path.isdir(path)):
            print(path,"is not the file")
        else:
            flexfilelist.append(path)
    flexfilelist.sort(key=str2int)
    print("\n".join(flexfilelist))
    for flexfile in flexfilelist:
        fpd=pd.read_table(flexfile,sep="\s+",header=None,usecols=[0,5,6],names=["chr","depth","GC"],dtype={"chr":int,"depth":int,"GC":float})
        
    #select bins
#         fpdfiltered=fpd
        fpd=fpd.iloc[selectedidx]
        fpdfiltered=fpd[fpd.chr<23].reset_index(drop=True)
        
        
    #depth_scale
        meandepth=np.mean(fpdfiltered.loc[:,"depth"])
        fpdfiltered.loc[:,"depth_scale"]=[e/meandepth for e in fpdfiltered.loc[:,"depth"]]# this step with warning, however I don't know why
    #correct the depth by GC 
        GCord_idx=list(fpdfiltered.loc[:,"GC"].sort_values().index)
        reGCord_idx=np.argsort(GCord_idx)
#         print(*fpdfiltered["depth_scale"].reindex(GCord_idx).values,sep="\n",file=open("temp_orddepth_scale",'w'))
#         print(*GCord_idx,file=open("GCord",'w'))
#         print(*reGCord_idx,sep="\n",file=open("temp_reord",'w'))
#         print("start spline")
        result=spline(fpdfiltered.loc[:,"depth_scale"].index.values,fpdfiltered["depth_scale"].reindex(GCord_idx).values)
#         print("finish spline",fpdfiltered.index.values.max())
#         xint=np.linspace(fpdfiltered.index.values.min(),fpdfiltered.index.values.max(),1000)
        print(*result(fpdfiltered.index.values),sep="\n",file=open("temp_spline.txt","w"))
        fpdfiltered["depth_correct"]=fpdfiltered["depth_scale"].values-result(fpdfiltered.index.values)[reGCord_idx]
        
        flexdatalist.append(fpdfiltered)

    for flex_idx1 in range(len(flexdatalist)):
        flex_data1=flexdatalist[flex_idx1]
#         print("\n",flexfilelist[flex_idx1],np.mean(d1["depth"]))
        for flex_idx2 in range(flex_idx1+1,len(flexdatalist)):
            flex_data2=flexdatalist[flex_idx2]
            print(flex_idx1,flex_idx2,end="\t")

            
    #calculate distance & file distance array
            a1=flex_data1.loc[:,"depth_scale"].values;b1=flex_data2.loc[:,"depth_scale"].values
            a2=flex_data1.loc[:,"depth_correct"].values;b2=flex_data2.loc[:,"depth_correct"].values
            
            if options.distype.lower()=="cos" or options.distype.lower()=="c":
                dist1=1-np.dot(a1,b1)/(np.linalg.norm(a1)*np.linalg.norm(b1))
                dist2=1-np.dot(a2,b2)/(np.linalg.norm(a2)*np.linalg.norm(b2))
            elif options.distype.lower()=="euclideanDistances" or options.distype.lower()=="e":
                dist1=np.linalg.norm(a1-b1)
                dist2=np.linalg.norm(a2-b2)
            elif options.distype.lower()=="hellinger" or options.distype.lower()=="h":
                dist1=1/np.sqrt(2)*np.linalg.norm(np.sqrt(a1)-np.sqrt(b1))
                dist2=1/np.sqrt(2)*np.linalg.norm(np.sqrt(a2)-np.sqrt(b2))
                print(a1,b1,dist1,a2,b2,dist2)
            mds_scale[flex_idx2][flex_idx1]=mds_scale[flex_idx1][flex_idx2]=str(dist1)
            mds_correct[flex_idx2][flex_idx1]=mds_correct[flex_idx1][flex_idx2]=str(dist2)
#     print("\nmds",mds)
    
    print("",end="\t",file=ofs)
    for fname in flexfilelist:
        print(re.search(r"[^/]*$",fname).group(0),end="\t",file=ofs)
    print("",file=ofs)
    n_i=0
    for distlist in mds_scale:
        print(re.search(r"[^/]*$",flexfilelist[n_i]).group(0),"\t".join(distlist),sep="\t",file=ofs)
        n_i+=1
    ofs.close()
    print("",end="\t",file=ofc)
    for fname in flexfilelist:
        print(re.search(r"[^/]*$",fname).group(0),end="\t",file=ofc)
    print("",file=ofc)
    n_i=0
    for distlist in mds_correct:
        print(re.search(r"[^/]*$",flexfilelist[n_i]).group(0),"\t".join(distlist),sep="\t",file=ofc)
        n_i+=1
    ofc.close()
    #