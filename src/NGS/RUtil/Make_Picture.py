import re
import rpy2.robjects as robjects

'''
Created on 2013-8-10

@author: rui
'''

class MakeMhtGraph(object):
    '''
    classdocs
    '''

    
    def __init__(self):
        super().__init__()
        self.originaldata={}
        self.dataForGraphe={}
        self.pathtoOutFileName=""
        '''
        Constructor
        '''
    def prepareMhtFile(self,inputfileName,dataType,chromPrefix="",postive_negtive=None):
        originalfile=open(inputfileName,'r')
        if postive_negtive==None:
            self.pathtoOutFileName=inputfileName+".z"+dataType
        else:
            self.pathtoOutFileName=inputfileName+"_"+postive_negtive+".z"+dataType
        for line in originalfile:
            linelist=re.split(r'\s+',line)
            currentChrom=linelist[0].strip()
            ChromNo=re.search(r"(\d+)$",currentChrom).group(1)
            if re.search(r"^"+chromPrefix,currentChrom):
                if linelist[4].strip()!="NA" or linelist[5].strip()!="NA":
                    if postive_negtive ==None:
                        if ChromNo in self.dataForGraphe.keys():
                            self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                        else:
                            self.dataForGraphe[ChromNo]=[tuple(linelist[1:])]
                    elif postive_negtive=="postive":
                        if float(linelist[5].strip())>0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                            else:
                                self.dataForGraphe[ChromNo]=[tuple(linelist[1:])]
                    elif postive_negtive=="negtive":
                        if float(linelist[5].strip())<0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                            else:
                                self.dataForGraphe[ChromNo]=[tuple(linelist[1:])]
                    else:
                        return "error"
        print(chromPrefix,"winNo","bp_start","bp_end",dataType,"z"+dataType,sep="\t",file=open(self.pathtoOutFileName,"w"))
#         print(chromPrefix,"bp_start",dataType)
        outfile=open(self.pathtoOutFileName,'a')
        for chromNo in sorted(self.dataForGraphe.keys(),key=lambda t:int(t)):
            for i in range(len(self.dataForGraphe[chromNo])):
                print(chromNo,*self.dataForGraphe[chromNo][i],sep="\t",file=outfile)
        outfile.close()
        return re.search(r"[^\\]*$", self.pathtoOutFileName).group(0)#for windows
    def makeMhtPictureFile(self,inputfileName,dataType,chromPrefix="",postive_negtive=None):
        name=self.prepareMhtFile(inputfileName,dataType,chromPrefix,postive_negtive)
        dir=re.search(r"^.*\\",self.pathtoOutFileName).group(0)
        r=robjects.r
        print(name,self.pathtoOutFileName,dir)
        
        r("setwd('"+dir+"')")
        r("library(gap)")
        
        r("tiff('"+name+".tiff')")
        
        r('x=read.table("'+self.pathtoOutFileName+'",header=T)')
        r("data=with(x,cbind("+chromPrefix+",bp_start,z"+dataType+"))")
        r('colors <- rep(c("green2","firebrick1"),38000)')
        r('par(las=1, xpd=TRUE, cex.axis=1.8, cex=0.5)')
        r('mhtplot(data,control=mht.control(logscale=FALSE,colors=colors,cex=4),pch=20,ylab="z'+dataType+'")')
        r('axis(2)')
        r('dev.off()')