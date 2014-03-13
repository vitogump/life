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
        self.originaldata = {}
        self.dataForGraphe = {}
        self.pathtoOutFileName = ""
        '''
        Constructor
        '''
    def prepareMhtFile(self, inputfileName, dataType, chromPrefix="", postive_negtive=None,fillvalue=0):
        originalfile = open(inputfileName, 'r')
        if postive_negtive == None:
            self.pathtoOutFileName = inputfileName + ".z" + dataType
        else:
            self.pathtoOutFileName = inputfileName + "_" + postive_negtive + ".z" + dataType
        for line in originalfile:
            linelist = re.split(r'\s+', line)
            currentChrom = linelist[0].strip()
            ChromNo = re.search(r"(\d+)$", currentChrom).group(1)
            if re.search(r"^" + chromPrefix, currentChrom):
                if linelist[4].strip() != "NA" or linelist[5].strip() != "NA" or True:
                    if postive_negtive == None:
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() != "NA":
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4]+[fillvalue,fillvalue]))
                        else:
                            self.dataForGraphe[ChromNo] = [tuple(linelist[1:4]+[fillvalue,fillvalue])]
#                             self.dataForGraphe[ChromNo] = [tuple(linelist[1:])]
                    elif postive_negtive == "postive":
                        if linelist[5].strip()=='NA' or float(linelist[5].strip()) <= 0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4]+[fillvalue,fillvalue]))
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4]+[fillvalue,fillvalue])]
                        else :#float(linelist[5].strip()) > 0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:])]
                    elif postive_negtive == "negtive":
                        if linelist[5].strip()=='NA' or float(linelist[5].strip()) >= 0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4]+[fillvalue,fillvalue]))
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4]+[fillvalue,fillvalue])]   
                        else:#float(linelist[5].strip()) < 0:
                            if ChromNo in self.dataForGraphe.keys():
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:]))
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:])]                            
                         
                    else:
                        return "error"
        print(chromPrefix, "winNo", "bp_start", "bp_end", dataType, "z" + dataType, sep="\t", file=open(self.pathtoOutFileName, "w"))
#         print(chromPrefix,"bp_start",dataType)
        outfile = open(self.pathtoOutFileName, 'a')
        for chromNo in sorted(self.dataForGraphe.keys(), key=lambda t:int(t)):
            for i in range(len(self.dataForGraphe[chromNo])):
                print(chromNo, *self.dataForGraphe[chromNo][i], sep="\t", file=outfile)
        outfile.close()
        return re.search(r"[^/]*$", self.pathtoOutFileName).group(0)  # for linux
    def makeMhtPicture_HistonPicture(self, inputfileName, dataType, chromPrefix="", postive_negtive=None,fillvalue=0):
        name = self.prepareMhtFile(inputfileName, dataType, chromPrefix, postive_negtive,fillvalue)
        dir = re.search(r"^.*/", self.pathtoOutFileName).group(0)
        r = robjects.r
        print(name, self.pathtoOutFileName, dir)
        
        r("setwd('" + dir + "')")
        r("library(gap)")
         
        r("pdf('" + name + ".pdf',width=10,height=4.5)")
         
        r('x=read.table("' + self.pathtoOutFileName + '",header=T)')
        r("data=with(x,cbind(" + chromPrefix + ",bp_start,z" + dataType + "))")
        r('colors <- rep(c("green2","firebrick1"),38000)')
        r('par(las=1, xpd=TRUE, cex.axis=1.0, cex=0.5)')
        r('mhtplot(data,control=mht.control(logscale=FALSE,colors=colors,cex=0.5),pch=20,ylab="z' + dataType + '")')
        r('axis(2)')
        r("title('"+name+"')")
        r('dev.off()')
#         print(name,dataType,"tiff('" + name + "histon_" + dataType + ".tiff'")
        r("tiff('" + name + "histon_" + dataType + ".tiff')")
#         if dataType == "Hp":
#             
#             r("bins=seq(0,0.6,by=0.001)")
#         elif dataType =="Fst":
#             r("bins=seq(-6,6,by=0.001)")
        r("hist(x$" + dataType + ",breaks=1000,main='"+name+"')")
        r('dev.off()')
        
        r("tiff('" + name + "histon_z" + dataType + ".tiff')")
#         if dataType == "Hp":
#             r("zbins=seq(-6,3.5,by=0.02)")
#         elif dataType == "Fst":
#             r("zbins=seq(-3.5,6,by=0.02)")
        
        r("hist(x$z" + dataType + ",breaks=1000,main='"+name+"')")
        r('dev.off()')

        
