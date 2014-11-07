'''
Created on 2013-8-10

@author: rui
'''

import itertools
import re, os, math

import numpy

from NGS.BasicUtil import Caculators
import NGS.BasicUtil.DerivedalleleProcessor as DAP
import rpy2.robjects as robjects


class Dstistics_allpop(object):
    def __init__(self, allpopslist):
        super().__init__()
        self.allpossiblecombination = list(itertools.permutations(allpopslist, 3))
    def caculateDofAllpossibleCombination(self,database,ip,usrname,pw,allpopssnptable, chromstable, winwidth, minlengthOfchrom, filenamepre):
        listtofinalfile = []
        D_sum_file=open(filenamepre+"D_SUM.txt","w")
        print("p1,p2,p3\tABBA\tBABA\tD-fixed\tSE-fixed\tD-SNP\tSE-SNP",file=D_sum_file)
        tempfiletomakebox = open(filenamepre + "test.box", 'w')
        
        print("D","group","chrom",sep='\t',file=tempfiletomakebox)
        for p1name, p2name, p3name in self.allpossiblecombination:
            allABBAcount = 0;allBABAcount = 0
            
            D = DAP.Dstatistics(database=database,ip=ip,usrname=usrname,pw=pw,allpopssnptable=allpopssnptable)
            D_caculator = Caculators.Caculate_Dstatistics()
            D.caculateFstAccordingdb(chromstable, p1name, p2name, p3name, D_caculator, winwidth, minlengthOfchrom)
            winCrossGenome_fix = []
            winCrossGenome_snp = []
            for chrom in sorted(D.DMapByChrom.keys()):
                if D.DMapByChrom[chrom][0][3] != 'NA':
                    winCrossGenome_snp.append(D.DMapByChrom[chrom][0][3])
                if D.DMapByChrom[chrom][0][2] != 'NA':
                    winCrossGenome_fix.append(D.DMapByChrom[chrom][0][2])
                allABBAcount += D.DMapByChrom[chrom][0][0]
                allBABAcount += D.DMapByChrom[chrom][0][1]
                print(str(D.DMapByChrom[chrom][0][3]),p1name+","+p2name+","+p3name,chrom,sep='\t',file=tempfiletomakebox)
            exception_fix = numpy.mean(winCrossGenome_fix);exception_snp = numpy.mean(winCrossGenome_snp)
            variance_fix = numpy.var(winCrossGenome_fix);variance_snp = numpy.var(winCrossGenome_snp)
            stderr_fix = math.sqrt(variance_fix * len(winCrossGenome_fix));stderr_snp = math.sqrt(variance_snp * len(winCrossGenome_snp))
            print(p1name +","+ p2name+"," + p3name,str(allABBAcount), str(allBABAcount), str(exception_fix), str(stderr_fix), str(exception_snp), str(stderr_snp),sep="\t",file=D_sum_file)
            listtofinalfile.append((p1name +","+ p2name+"," + p3name, str(allABBAcount), str(allBABAcount), str(exception_fix), str(stderr_fix), str(exception_snp), str(stderr_snp)))
            D.dbtools.disconnect()
        tempfiletomakebox.close()
#         for rec in listtofinalfile:
#             print("\t".join(rec),file=D_sum_file)
        D_sum_file.close()
        

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
    def prepareMhtFileWithgeneName(self, inputfileName, dataType, chromPrefix="", positive_negtive=None, fillvalue=0):
        """fill all NA value window with fillvalue,and fill all window that zvalue<=0 with fillvalue when positive_negtive= positive....
        """
        originalfile = open(inputfileName, 'r')
        
        print("title", originalfile.readline())
        if positive_negtive == None:
            print(inputfileName, dataType)
            self.pathtoOutFileName = inputfileName + ".z" + dataType
        else:
            self.pathtoOutFileName = inputfileName + "_" + positive_negtive + ".z" + dataType
        for line in originalfile:
            linelist = re.split(r'\s+', line.strip())
            currentChrom = linelist[0].strip()
            try:
                ChromNo = re.search(r"([\d]+)", currentChrom).group(1)
            except AttributeError:
                ChromString = re.search(r"" + chromPrefix + "([\d\D]+)$", currentChrom).group(1)
                ChromNo = 0
                for e in ChromString:
                    if e.isalpha():
                        ChromNo += ord(e)
                    elif e.isdigit():
                        ChromNo += int(e) 
                print(ChromNo)
            if re.search(r"^" + chromPrefix, currentChrom):
                if linelist[4].strip() != "NA" or linelist[5].strip() != "NA" or True:
                    if positive_negtive == None:
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() != "NA":
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:7]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]]))
                        else:
                            self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]])]
                    elif positive_negtive == "positive":
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) <= 0:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:7]))
                        else:
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) <= 0:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]])]
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:7])]
                    elif positive_negtive == "negtive":
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) >= 0:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:7]))
                        else:
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) >= 0:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]])]
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:7])]                           
                    else:
                        return "error"
                else:
                    if ChromNo in self.dataForGraphe.keys():
                        self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue, linelist[6]]))
                    
#         if len(linelist)==6:
#             print(chromPrefix, "winNo", "bp_start", "bp_end", dataType, "z" + dataType, sep="\t", file=open(self.pathtoOutFileName, "w"))
        print(chromPrefix, "winNo", "bp_start", "bp_end", dataType, "z" + dataType, "geneName", sep="\t", file=open(self.pathtoOutFileName, "w"))
#         print(chromPrefix,"bp_start",dataType)
        outfile = open(self.pathtoOutFileName, 'a')
        for chromNo in sorted(self.dataForGraphe.keys(), key=lambda t:int(t)):
            for i in range(len(self.dataForGraphe[chromNo])):
                print(chromNo, *self.dataForGraphe[chromNo][i], sep="\t", file=outfile)
        outfile.close()
        return re.search(r"[^/]*$", self.pathtoOutFileName).group(0)  # for linux


    def prepareMhtFile(self, inputfileName, dataType, chromPrefix="", positive_negtive=None, fillvalue=0):
        """fill all NA value window with fillvalue,and fill all window that zvalue<=0 with fillvalue when positive_negtive= positive....
        """
        originalfile = open(inputfileName, 'r')
        print("title", originalfile.readline())
        if positive_negtive == None:
            print(inputfileName, dataType)
            self.pathtoOutFileName = inputfileName + ".z" + dataType
        else:
            self.pathtoOutFileName = inputfileName + "_" + positive_negtive + ".z" + dataType
        for line in originalfile:
            linelist = re.split(r'\s+', line.strip())
            currentChrom = linelist[0].strip()
            try:
                ChromNo = re.search(r"([\d]+)", currentChrom).group(1)
            except AttributeError:
                ChromString = re.search(r"" + chromPrefix + "([\d\D]+)$", currentChrom).group(1)
                ChromNo = 0
                for e in ChromString:
                    if e.isalpha():
                        ChromNo += ord(e)
                    elif e.isdigit():
                        ChromNo += int(e) 
                print(ChromNo)
            if re.search(r"^" + chromPrefix, currentChrom):
                if linelist[4].strip() != "NA" or linelist[5].strip() != "NA" or True:
                    if positive_negtive == None:
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() != "NA":
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:6]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue]))
                        else:
                            self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue])]
                    elif positive_negtive == "positive":
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) <= 0:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:6]))
                        else:
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) <= 0:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue])]
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:6])]
                    elif positive_negtive == "negtive":
                        if ChromNo in self.dataForGraphe.keys():
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) >= 0:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue]))
                            else:
                                self.dataForGraphe[ChromNo].append(tuple(linelist[1:6]))
                        else:
                            if linelist[5].strip() == 'NA' or float(linelist[5].strip()) >= 0:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:4] + [fillvalue, fillvalue])]
                            else:
                                self.dataForGraphe[ChromNo] = [tuple(linelist[1:6])]                           
                    else:
                        return "error"
                else:
                    if ChromNo in self.dataForGraphe.keys():
                        self.dataForGraphe[ChromNo].append(tuple(linelist[1:4] + [fillvalue, fillvalue]))
                    
#         if len(linelist)==6:
#             print(chromPrefix, "winNo", "bp_start", "bp_end", dataType, "z" + dataType, sep="\t", file=open(self.pathtoOutFileName, "w"))
        print(chromPrefix, "winNo", "bp_start", "bp_end", dataType, "z" + dataType, sep="\t", file=open(self.pathtoOutFileName, "w"))
#         print(chromPrefix,"bp_start",dataType)
        outfile = open(self.pathtoOutFileName, 'a')
        for chromNo in sorted(self.dataForGraphe.keys(), key=lambda t:int(t)):
            for i in range(len(self.dataForGraphe[chromNo])):
                print(chromNo, *self.dataForGraphe[chromNo][i], sep="\t", file=outfile)
        outfile.close()
        return re.search(r"[^/]*$", self.pathtoOutFileName).group(0)  # for linux
    def makeMhtPicture_HistonPicture(self, inputfileName, dataType, chromPrefix="", positive_negtive=None, fillvalue=0):
        line = open(inputfileName, 'r').readline().strip()
        open(inputfileName, 'r').close()
        if 8 == len(re.split(r'\s+', line)):
            name = self.prepareMhtFileWithgeneName(inputfileName, dataType, chromPrefix, positive_negtive, fillvalue)
        elif 6 == len(re.split(r'\s+', line)):
            name = self.prepareMhtFile(inputfileName, dataType, chromPrefix, positive_negtive, fillvalue)
        dir = re.search(r"^.*/", self.pathtoOutFileName).group(0)
        r = robjects.r
        print(name, self.pathtoOutFileName, dir)
        
        r("setwd('" + dir + "')")
        r("library(gap)")
        r("library(Cairo)")
         
        r("pdf('" + name + ".pdf',width=10,height=4.5)")
        print('x=read.table("' + self.pathtoOutFileName + '",header=T)') 
        r('x=read.table("' + self.pathtoOutFileName + '",header=T)')
        
        r("data=with(x,cbind(" + chromPrefix + ",bp_start,z" + dataType + "))")
        r('colors <- rep(c("green2","firebrick1"),38000)')
        r('par(las=1, xpd=TRUE, cex.axis=1.0, cex=0.5)')
        r('mhtplot(data,control=mht.control(logscale=FALSE,colors=colors,cex=0.5),pch=20,ylab="z' + dataType + '")')
        r('axis(2)')
        r("title('" + name + "')")
        r('dev.off()')
#         print(name,dataType,"tiff('" + name + "histon_" + dataType + ".tiff'")
        r("CairoPNG('" + name + "histon_" + dataType + ".png')")
#         if dataType == "Hp":
#             
#             r("bins=seq(0,0.6,by=0.001)")
#         elif dataType =="Fst":
#             r("bins=seq(-6,6,by=0.001)")
        r("hist(x$" + dataType + ",breaks=1000,main='" + name + "')")
        r('dev.off()')
        
        r("CairoPNG('" + name + "histon_z" + dataType + ".png')")
#         if dataType == "Hp":
#             r("zbins=seq(-6,3.5,by=0.02)")
#         elif dataType == "Fst":
#             r("zbins=seq(-3.5,6,by=0.02)")
        
        r("hist(x$z" + dataType + ",breaks=1000,main='" + name + "')")
        r('dev.off()')

        
