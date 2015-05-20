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
# import rpy2.rinterface as ri


# ri.set_initoptions((b'rpy2', b'--verbose', b'--no-save'))
# ri.initr()

class Dstistics_allpop(object):
    def __init__(self, allpopslist):
        super().__init__()
        self.allpossiblecombination = list(itertools.permutations(allpopslist, 3))
    def caculateDofAllpossibleCombination(self, database, ip, usrname, pw, allpopssnptable, chromstable, winwidth, minlengthOfchrom, filenamepre):
        listtofinalfile = []
        D_sum_file = open(filenamepre + "D_SUM.txt", "w")
        print("p1,p2,p3\tABBA\tBABA\tD-fixed\tSE-fixed\tD-SNP\tSE-SNP", file=D_sum_file)
        tempfiletomakebox = open(filenamepre + "test.box", 'w')
        
        print("D", "group", "chrom", sep='\t', file=tempfiletomakebox)
        for p1name, p2name, p3name in self.allpossiblecombination:
            allABBAcount = 0;allBABAcount = 0;noofsnp = 0
            
            D = DAP.Dstatistics(database=database, ip=ip, usrname=usrname, pw=pw, allpopssnptable=allpopssnptable)
            D_caculator = Caculators.Caculate_Dstatistics()
            D.caculateFstAccordingdb(chromstable, p1name, p2name, p3name, D_caculator, winwidth, minlengthOfchrom)
            winCrossGenome_fix = []
            winCrossGenome_snp = []
            for chrom in sorted(D.DMapByChrom.keys()):
                if D.DMapByChrom[chrom][0][4] != 'NA':
                    winCrossGenome_snp.append(D.DMapByChrom[chrom][0][4])
                if D.DMapByChrom[chrom][0][2] != 'NA':
                    winCrossGenome_fix.append(D.DMapByChrom[chrom][0][2])
                allABBAcount += D.DMapByChrom[chrom][0][0]
                allBABAcount += D.DMapByChrom[chrom][0][1]
                noofsnp += D.DMapByChrom[chrom][0][3]
                print(str(D.DMapByChrom[chrom][0][3]), str(D.DMapByChrom[chrom][0][4]), p1name + "," + p2name + "," + p3name, chrom, sep='\t', file=tempfiletomakebox)
            exception_fix = numpy.mean(winCrossGenome_fix);exception_snp = numpy.mean(winCrossGenome_snp)
            variance_fix = numpy.var(winCrossGenome_fix);variance_snp = numpy.var(winCrossGenome_snp)
            stderr_fix = math.sqrt(variance_fix * len(winCrossGenome_fix));stderr_snp = math.sqrt(variance_snp * len(winCrossGenome_snp))
            print(p1name + "," + p2name + "," + p3name, str(allABBAcount), str(allBABAcount), str(exception_fix), str(stderr_fix), str(noofsnp), str(exception_snp), str(stderr_snp), sep="\t", file=D_sum_file)
            listtofinalfile.append((p1name + "," + p2name + "," + p3name, str(allABBAcount), str(allBABAcount), str(exception_fix), str(stderr_fix), str(noofsnp), str(exception_snp), str(stderr_snp)))
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
        a = os.popen("pwd")
        self.olddir = a.readline().strip()
        a.close()
        self.originaldata = {}
        self.dataForGraphe = {}
        self.pathtoOutFileName = ""
        self.namewithoutpath=""
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
    def prepareMhtFile(self,inputfileName, dataType, positive_negtive=None, fillvalue=0):
        pathtoOutFileName_new = inputfileName + "_" + positive_negtive + ".z" + dataType
        fillvalue = str(fillvalue)
        if not fillvalue.isnumeric():
            fillvalue='"'+fillvalue+'"'
        if positive_negtive == "positive":
            os.system(""" awk '{OFS="\\t";if(NR!=1 && ($6=="NA" || $7<0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
        elif positive_negtive == "negtive":
            print(""" awk '{OFS="\\t";if(NR!=1 && ($6=="NA" || $7>0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
            os.system(""" awk '{OFS="\\t";if(NR!=1 && ($6=="NA" || $7>0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
        return re.search(r"[^/]*$", pathtoOutFileName_new).group(0),pathtoOutFileName_new  # for linux
    
    def makeHistonPicture(self,inputfileName, dataType,ylim=None,xlim=None,columnnames=("winvalue","zvalue")):
        r = robjects.r
        if re.search(r"^.*/", inputfileName) != None:
            dir = re.search(r"^.*/", inputfileName).group(0)
        else:
            dir = self.olddir
        
        namewithoutpath=re.search(r"[^/]*$", inputfileName).group(0)
        r("setwd('" + dir + "')")
        r('.libPaths("/opt/Rpackages/")')
        r("library(Cairo)")
        r('x=read.table("' + namewithoutpath + '",header=T)')
        print(dir,namewithoutpath)
        if ylim!=None and xlim!=None:
            print("=========================================================================================================================================")
            for columnname in columnnames:
                r("CairoPNG('" + namewithoutpath + "histon_"+columnname+"_" + dataType + ".png',width=1600,height=800)")
                r("hist(x$"+columnname+",breaks=500,ylim=c(0,6000),xlim=c(0,60),main='"+ namewithoutpath +"',ylab='No.bins',xlab='SNPs per Kb')")
                sd=r("sd(x$"+columnname+",na.rm=TRUE)")
                mean=r("mean(x$"+columnname+",na.rm=TRUE)")
                print(inputfileName,columnname,sd,mean,file=open(dir+"test.txt",'a'))
                r('dev.off()')
        else:
            print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            for columnname in columnnames:
                r("CairoPNG('" + namewithoutpath + "histon_"+columnname+"_" + dataType + ".png',width=1600,height=800)")
                r("hist(x$"+columnname+",breaks=1000,main='" + namewithoutpath + "')")
                sd=r("sd(x$"+columnname+",na.rm=TRUE)")
                mean=r("mean(x$"+columnname+",na.rm=TRUE)")
#                 print(inputfileName,columnname,sd,mean,file=open(dir+"/test.txt",'a'))
                r('dev.off()')
        os.system("cd "+self.olddir)
        print(ylim,xlim)     
    def makeMhtplots_compareInOnePicture_withgeneName(self, outputnamewithpath,positive_winfiles,negtive_winfiles,fillvalue=0,columnname="zvalue"):
        scriptfile=open("stripts.R",'w')
        print(outputnamewithpath,file=scriptfile)
        print(positive_winfiles,negtive_winfiles)
        if re.search(r"^.*/", outputnamewithpath)!=None:
            dir = re.search(r"^.*/", outputnamewithpath).group(0)
        else:
            dir=self.olddir
        outname=re.search(r"[^/]*$", outputnamewithpath).group(0)
        r = robjects.r
        positive_filenames=[""]*len(positive_winfiles);positive_filenameWithPaths=[""]*len(positive_winfiles)
        negtive_filenames=[""]*len(negtive_winfiles);negtive_filenameWithPaths=[""]*len(negtive_winfiles)
        p_threshold=[0]*len(positive_winfiles)
        n_threshold=[0]*len(negtive_winfiles)
        r("setwd('" + dir + "')")
        r('.libPaths("/opt/Rpackages/")')
        r("library(gap)")
        r("library(Cairo)")
#         r('CairoPNG("'+outname+'.png",width='+str(((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)*2)+',height='+str((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)+')')
        r('CairoPNG("'+outname+'.png",width=1600,height=800)')
        for i in range(0,len(positive_winfiles)):
            p_threshold[i]=float(positive_winfiles[i][1].strip())
            positive_filenames[i],positive_filenameWithPaths[i]=self.prepareMhtFile(positive_winfiles[i][0], "Fst", "positive", fillvalue)
            print('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)',file=scriptfile)
            r('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)')
            print('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
            r('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
            print('p_highlight'+str(i)+'<- p_dataframe'+str(i)+'[p_dataframe'+str(i)+'$'+columnname+'>='+str(p_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]',file=scriptfile)
            r('p_highlight'+str(i)+'<- p_dataframe'+str(i)+'[p_dataframe'+str(i)+'$'+columnname+'>='+str(p_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]')
            print('p_highlithcolors'+str(i)+'<- rep("red",nrow(p_highlight'+str(i)+'))',file=scriptfile)
            r('p_highlithcolors'+str(i)+'<- rep("red",nrow(p_highlight'+str(i)+'))')
        for i in range(0,len(negtive_winfiles)):
            n_threshold[i]=float(negtive_winfiles[i][1].strip())
            negtive_filenames[i],negtive_filenameWithPaths[i]=self.prepareMhtFile(negtive_winfiles[i][0], "Hp", "negtive", fillvalue)
            r('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)')
            r('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
            print('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)',file=scriptfile)
            print('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
            r('n_highlight'+str(i)+'<-n_dataframe'+str(i)+'[n_dataframe'+str(i)+'$'+columnname+'<='+str(n_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]')
            print('n_highlight'+str(i)+'<-n_dataframe'+str(i)+'[n_dataframe'+str(i)+'$'+columnname+'<='+str(n_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]',file=scriptfile)
            print('n_highlithcolors'+str(i)+'<- rep("red",nrow(n_highlight'+str(i)+'))',file=scriptfile)
            r('n_highlithcolors'+str(i)+'<- rep("red",nrow(n_highlight'+str(i)+'))')
            
        r('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red"),300)')
        r('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(0.8, 4, 0.8, 2))')
        print('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red"),300)',file=scriptfile)
        print('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(0.8, 4, 0.8, 2))',file=scriptfile)
        if len(positive_winfiles)==1:
            hopscex='1'
        else:
            hopscex="0.65"
        for i in range(0,len(positive_winfiles)):
            r('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)')
            
            r('hops<-hmht.control(data=p_highlight'+str(i)+',cex='+hopscex+')')
            r('mhtplot(p_data'+str(i)+',ops,hops,pch=19,ylab="z' + "Fst" + '",xlab="")')
            r("title(main='" + positive_filenames[i] + "',cex.main=0.8)")
            r('axis(2)')
            r('abline(h=5)')
            print('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)',file=scriptfile)
            print('hops<-hmht.control(data=p_highlight'+str(i)+',cex=0.6)',file=scriptfile)
            print('mhtplot(p_data'+str(i)+',ops,hops,pch=19,ylab="z' + "Fst" + '",xlab="")',file=scriptfile)
            print("title(main='" + positive_filenames[i] + "',cex.main=0.8)",file=scriptfile)
            print('axis(2)',file=scriptfile)
            print('abline(h=5)',file=scriptfile)
        for i in range(0,len(negtive_winfiles)):
            r('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)')
            r('hops<-hmht.control(data=n_highlight'+str(i)+',cex=0.6)')
            r('mhtplot(n_data'+str(i)+',ops,hops,pch=19,ylab="zHp",xlab="")')
            r("title(main='" + negtive_filenames[i] + "',cex.main=0.8)")
            r('axis(2)')
            r('abline(h=-5)')
            print('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)',file=scriptfile)
            print('hops<-hmht.control(data=n_highlight'+str(i)+',cex=0.6)',file=scriptfile)
            print('mhtplot(n_data'+str(i)+',ops,hops,pch=19,ylab="zHp",xlab="")',file=scriptfile)
            print("title(main='" + negtive_filenames[i] + "',cex.main=0.8)",file=scriptfile)
            print('axis(2)',file=scriptfile)
            print('abline(h=-5)',file=scriptfile)
        r('axis(1)')
        r('dev.off()')
        print(r('Cairo.capabilities()'))
        scriptfile.close()

    def makeMhtplots_compareInOnePicture(self, outputnamewithpath,positive_winfiles,negtive_winfiles,fillvalue=0,columnname="zvalue"):
        scriptfile=open("stripts.R",'w')
        print(outputnamewithpath,file=scriptfile)
        print(positive_winfiles,negtive_winfiles)
        if re.search(r"^.*/", outputnamewithpath)!=None:
            dir = re.search(r"^.*/", outputnamewithpath).group(0)
        else:
            dir=self.olddir
        outname=re.search(r"[^/]*$", outputnamewithpath).group(0)
        r = robjects.r
        positive_filenames=[""]*len(positive_winfiles);positive_filenameWithPaths=[""]*len(positive_winfiles)
        negtive_filenames=[""]*len(negtive_winfiles);negtive_filenameWithPaths=[""]*len(negtive_winfiles)
        r("setwd('" + dir + "')")
        r('.libPaths("/opt/Rpackages/")')
        r("library(gap)")
        r("library(Cairo)")
#         r('CairoPNG("'+outname+'.png",width='+str(((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)*2)+',height='+str((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)+')')
        r('CairoPNG("'+outname+'.png",width=1600,height=800)')
        for i in range(0,len(positive_winfiles)):
            positive_filenames[i],positive_filenameWithPaths[i]=self.prepareMhtFile(positive_winfiles[i], "Fst", "positive", fillvalue)
            print('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
            r('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)')
            r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))')
            print('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
#             r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,zvalue))')
        for i in range(0,len(negtive_winfiles)):
            negtive_filenames[i],negtive_filenameWithPaths[i]=self.prepareMhtFile(negtive_winfiles[i], "Hp", "negtive", fillvalue)
            r('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)')
            r('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
            'n_data'+str(i)+' <- with(n_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))'
            print('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
            print('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
        r('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red"),300)')
        r('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(2, 4, 1.5, 2))')
        print('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red"),300)',file=scriptfile)
        print('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(2, 4, 1.5, 2))',file=scriptfile)
        for i in range(0,len(positive_winfiles)):
            r('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=16,ylab="z' + "Fst" + '",xlab="")')
            r("title(main='" + positive_filenames[i] + "',cex.main=2)")
            r('axis(2)')
            print('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=16,ylab="z' + "Fst" + '",xlab="")',file=scriptfile)
            print("title(main='" + positive_filenames[i] + "',cex.main=2)",file=scriptfile)
            print('axis(2)',file=scriptfile)
        for i in range(0,len(negtive_winfiles)):
            r('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=16,ylab="zHp",xlab="")')
            r("title(main='" + negtive_filenames[i] + "',cex.main=2)")
            r('axis(2)')
            print('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=16,ylab="zHp",xlab="")',file=scriptfile)
            print("title(main='" + negtive_filenames[i] + "',cex.main=2)",file=scriptfile)
            print('axis(2)',file=scriptfile)
        r('axis(1)')
        r('dev.off()')
        print(r('Cairo.capabilities()'))
        scriptfile.close()
        
