# -*- coding: utf-8 -*- 
'''
Created on 2014-11-8

@author: liurui
'''
import os,re

class OperatorWithData():
    def __init__(self):
        pass
    def process(self,p,d):
        print(p,d)
#myprint=OperatorWithData()
class OperatorWithData_mode1(OperatorWithData):
    def __init__(self,cmdline):
        super().__init__()
        self.cmdline=cmdline
    def process(self,curpath,curdepth):
        print("mode1 process")
        newcmdline=self.cmdline
        subtargets=re.findall(r"\${.*?}",newcmdline)
        targetdatasuffix=[]
        for target in subtargets:
            c=re.search(r'\${(.*?)}',target).group(1)
            if "output"==c:
                print(target,subtargets)
                subtargets.remove(target)
                continue
            targetdatasuffix.append(c)
        
        datafiles=os.listdir(path=curpath)
        print("OperatorWithData_mode1",datafiles)
        for i in range(0,len(targetdatasuffix)):
            for datafilename in datafiles:
                if re.search(r".*?"+targetdatasuffix[i],datafilename)!=None:
                    newcmdline=re.sub(r"\${\s*"+targetdatasuffix[i]+"\s*}"," "+curpath+"/"+datafilename+" ",newcmdline)
                    # sub was acted from the first to the rear most
        return newcmdline
class OperatorWithData_mode2(OperatorWithData):
    def __init__(self,cmdline):
        super().__init__()
        self.cmdline=cmdline
        self.newcmdline=""
    def process(self,curpath,curdepth):
        print("mode2 process")
        newcmdline=self.cmdline

        option_suffix_obj=re.search(r"(-\wd+)\s+\${(.*?)}",newcmdline)# for example "-INPUT ${.bam}"
        optionstr=option_suffix_obj.group(1)
        suffixstr=option_suffix_obj.group(2)
        
        datafiles=os.listdir(path=curpath)
        print("OperatorWithData_mode2",datafiles)
        for datafilename in datafiles:
            if re.search(r".*?"+suffixstr,datafilename)!=None:
                newcmdline=re.sub(r"-\wd+\s+\${.*?}",curpath+"/"+datafilename+" "+option_suffix_obj.group(0),newcmdline)

        self.newcmdline=newcmdline
        return newcmdline
def upTodownTravelDir(rootDir,OperatorWithData,maxdepth=9999,curdepth=0,Interceptor=None):
    """
        Interceptor=([subdir names list],depth of the names expected)
    """
    print(rootDir)
    if maxdepth==0:
        newcmdline=OperatorWithData.process(rootDir,curdepth)
        print(rootDir,newcmdline)
        return
    # now go into a deeper dir
    curdepth=curdepth+1
    maxdepth=maxdepth-1
    for elem in os.listdir(path=rootDir):
        path = rootDir+"/"+elem
        if not os.path.isdir(path):
            #this is a data file
            pass
            #print("data",path)
        else:
            #this is a folder
            if Interceptor!=None and curdepth==Interceptor[1] and elem not in Interceptor[0]:
                continue 
            print("go into folder",path)
            upTodownTravelDir(path,OperatorWithData,maxdepth,curdepth)