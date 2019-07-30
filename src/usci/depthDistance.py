# -*- coding: UTF-8 -*-
'''
Created on 2019��7��30��

@author: liurui
'''

from optparse import OptionParser
import os,re
import pandas as pd
import numpy as np

parser = OptionParser()
parser.add_option("-d", "--flexDir", dest="flexDir", help="")
parser.add_option("-s", "--selectbinfiles", dest="selectbinfiles")

parser.add_option("-o", "--output", dest="output")

(options, args) = parser.parse_args()
flexDir=options.flexDir
flexfilelist=[]
f=open(options.selectbinfiles,'r')
selectedidx=[int(e.strip()) for e in f.readlines()];f.close()

def tryint(s):
    try:
        print(int(s))
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
mds=[]
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

    for flex_idx1 in range(len(flexfilelist)):
        flex_data1=pd.read_table(flexfilelist[flex_idx1],sep="\s+",header=None,usecols=[0,5,6],names=["chr","depth","GC"],dtype={"chr":int,"depth":int,"GC":float})
        
#         print("\n",flexfilelist[flex_idx1],np.mean(d1["depth"]))
        mds.append([])
        for flex_idx2 in range(flex_idx1+1,len(flexfilelist)):
            flex_data2=pd.read_table(flexfilelist[flex_idx2],sep="\s+",header=None,usecols=[0,5,6],names=["chr","depth","GC"],dtype={"chr":int,"depth":int,"GC":float})
            
            print(flex_idx1,flex_idx2,end="\t")
    #select bins
            d1=flex_data1.iloc[selectedidx];d1=d1[d1.chr<23]
            d2=flex_data2.iloc[selectedidx];d2=d2[d2.chr<23]
            d1["depth_scale"]=[e/np.mean(d1["depth"]) for e in d1["depth"]]
            d2["depth_scale"]=[e/np.mean(d2["depth"]) for e in d2["depth"]]
    #correct the depth by GC 
            
    #calculate distance & file distance array
            a=d1["depth_scale"].as_matrix();b=d2["depth_scale"].as_matrix()
            dist=np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))
            mds[-1].append(str(dist))
    print(mds)
    for distlist in mds:
        print("\t".join(distlist),file=of)
    of.close()
    #