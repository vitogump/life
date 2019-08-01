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

parser = OptionParser()
parser.add_option("-d", "--flexDir", dest="flexDir", help="")
parser.add_option("-s", "--selectbinfiles", dest="selectbinfiles")

parser.add_option("-o", "--output", dest="output")

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

mds=[["0" for i in range(75)] for j in range(75)]
of=open(options.output,"w")
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
        fpd=fpd.iloc[selectedidx]
        fpd23=fpd[fpd.chr<23]
        
    #depth_scale
        meandepth=np.mean(fpd23.loc[:,"depth"])
        fpd23.loc[:,"depth_scale"]=[e/meandepth for e in fpd23.loc[:,"depth"]]# this step with warning, however I don't know why
    #correct the depth by GC 
        fpd23orderd=fpd23.loc[:,"GC"].sort_values()
        GCord_idx=list(fpd23orderd.index)
        result=spline(fpd23.loc[:,"depth"].index.values,fpd23orderd["depth"])
        result(fpd23.loc[:,"depth"].index.values)
        print(type(result))
        flexdatalist.append(fpd23)
    for flex_idx1 in range(len(flexdatalist)):
        flex_data1=flexdatalist[flex_idx1]
#         print("\n",flexfilelist[flex_idx1],np.mean(d1["depth"]))
        for flex_idx2 in range(flex_idx1+1,len(flexdatalist)):
            flex_data2=flexdatalist[flex_idx2]
            print(flex_idx1,flex_idx2,end="\t")

            
    #calculate distance & file distance array
            a=flex_data1.loc[:,"depth_scale"].values;b=flex_data2.loc[:,"depth_scale"].values
            dist=np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))
            mds[flex_idx2][flex_idx1]=mds[flex_idx1][flex_idx2]=str(dist)
            
#     print("\nmds",mds)
    
    print("",end="\t",file=of)
    for fname in flexfilelist:
        print(re.search(r"[^/]*$",fname).group(0),end="\t",file=of)
    print("",file=of)
    n_i=0
    for distlist in mds:
        print(re.search(r"[^/]*$",flexfilelist[n_i]).group(0),"\t".join(distlist),sep="\t",file=of)
        n_i+=1
    of.close()
    #