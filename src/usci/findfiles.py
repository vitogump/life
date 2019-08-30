'''
Created on 2019年8月26日

@author: liurui
'''
import re,os
import pandas as pd


if __name__ == '__main__':
    of=open("fqfileTrace.txt",'w')
    fpd=pd.read_excel('lung_breast_334_filter.xls',sheet_name='lung_breast_334_info')
    for s_idx in fpd.index.values:
        plasma_sample,plasma_flex=fpd.ix[s_idx,['plasma_sample','plasma_flex']].values
        flexpath=re.search(r"([\w\W]*)/[^/]*$",plasma_flex).group(1)
        for elem in os.listdir(flexpath+os.sep+"../result/sh"):
            if (not os.path.isdir(flexpath+os.sep+"../result/sh/"+elem)) and re.search(r""+plasma_sample,elem)!=None and "NIPT1" in elem:
                shfile=open(flexpath+os.sep+"../result/sh/"+elem,'r')
                for line in shfile:
                    if "R1.fq.gz" in line:
                        for pathtofq in re.split(r'\s+',line.strip()):
                            if "R1.fq.gz" in pathtofq:pathtofq=re.search(r"([\w\W]*)/[^/]*R1.fq.gz$",pathtofq).group(1).replace("data","project",1);break
                        break           
                print(plasma_sample,pathtofq,plasma_flex,re.search(r'new/([^/]*)',plasma_flex).group(1),end="\t",sep="\t",file=of)
                break
        else:
            print(plasma_sample,plasma_flex,"not find",file=of)
        
        leucocyte_sample,leucocyte_flex=fpd.ix[s_idx,['leucocyte_sample','leucocyte_flex']].values
#         print(leucocyte_sample,leucocyte_flex)
        if "/flex" in leucocyte_flex:flexpath=re.search(r"([\w\W]*)/[^/]*$",leucocyte_flex).group(1)
        for elem in os.listdir(flexpath+os.sep+"../result/sh"):
            if (not os.path.isdir(flexpath+os.sep+"../result/sh/"+elem)) and re.search(r""+leucocyte_sample,elem)!=None and "NIPT1" in elem:
                shfile=open(flexpath+os.sep+"../result/sh/"+elem,'r')
                for line in shfile:
                    if "R1.fq.gz" in line:
                        for pathtofq in re.split(r'\s+',line.strip()):
                            if "R1.fq.gz" in pathtofq:pathtofq=re.search(r"([\w\W]*)/[^/]*R1.fq.gz$",pathtofq).group(1).replace("data","project",1);break
                            
                        break
                print(leucocyte_sample,pathtofq,leucocyte_flex,re.search(r'new/([^/]*)',leucocyte_flex).group(1),sep="\t",file=of)
                break
        else:
            print(leucocyte_sample,leucocyte_flex,"not find",file=of)
    of.close()
#         print(pathtoflex,plasma_sample)
