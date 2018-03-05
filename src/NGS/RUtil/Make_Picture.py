'''
Created on 2013-8-10

@author: rui
'''

import itertools
from os.path import sys
import re, os, math, time

import numpy

from NGS.BasicUtil import Caculators, Util
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
            D.caculateDstatisticsAccordingdb(chromstable, p1name, p2name, p3name, D_caculator, winwidth, minlengthOfchrom)
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
            print(""" awk '{OFS="\\t";if(NR!=1 && ($7=="NA" || $7<0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
            os.system(""" awk '{OFS="\\t";if(NR!=1 && ($7=="NA" || $7<0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
        elif positive_negtive == "negtive":
            print(""" awk '{OFS="\\t";if(NR!=1 && ($7=="NA" || $7>0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
            os.system(""" awk '{OFS="\\t";if(NR!=1 && ($7=="NA" || $7>0)){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
        elif positive_negtive=="allvalue":
            os.system(""" awk '{OFS="\\t";if(NR!=1 && ($7=="NA" )){$6=""" + fillvalue + """;$7=""" + fillvalue + """};print $0}' """ + inputfileName + ">" + pathtoOutFileName_new)
        print(pathtoOutFileName_new)
        return re.search(r"[^/]*$", pathtoOutFileName_new).group(0),pathtoOutFileName_new  # for linux
    
    def makeHistonPicture(self,inputfileName, dataType,ylim=None,xlim=None,columnnames=("winvalue","zvalue")):
        r = robjects.r
        if re.search(r"^.*/", inputfileName) != None:
            dir = re.search(r"^.*/", inputfileName).group(0)
        else:
            dir = self.olddir
        
        namewithoutpath=re.search(r"[^/]*$", inputfileName).group(0)
        r("setwd('" + dir + "')")
        r('.libPaths("/home/liurui/software/Rpackages")')
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
                try:
                    r("hist(x$"+columnname+",breaks=1000,main='" + namewithoutpath + "')")
                    sd=r("sd(x$"+columnname+",na.rm=TRUE)")
                    mean=r("mean(x$"+columnname+",na.rm=TRUE)")
                except:
                    pass
#                 print(inputfileName,columnname,sd,mean,file=open(dir+"/test.txt",'a'))
                r('dev.off()')
        os.system("cd "+self.olddir)
        print(ylim,xlim,"mkmhthist done")     
    def makeMhtplots_compareInOnePicture_withgeneName(self, outputnamewithpath,positive_winfiles,negtive_winfiles,fillvalue=0,columnname="zvalue"):
        scriptfile=open("stripts_withgene.R",'w')
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
        r('.libPaths("/home/liurui/software/Rpackages")')
        r("library(gap)")
        r("library(Cairo)")
#         r('CairoPNG("'+outname+'.png",width='+str(((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)*2)+',height='+str((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)+')')
        r('CairoPNG("'+outname+'.png",width=1800,height=800)')
#         r('CairoPDF("'+outname+'.pdf"')
        for i in range(0,len(positive_winfiles)):
            threshold_title_list=re.split(r"_",positive_winfiles[i][1].strip())
            p_threshold[i]=float(threshold_title_list[0])
            positive_filenames[i],positive_filenameWithPaths[i]=self.prepareMhtFile(positive_winfiles[i][0], "Fst", "positive", fillvalue)
            print('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)',file=scriptfile)
            r('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)')
            print('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
            r('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
            print('p_highlight'+str(i)+'<- p_dataframe'+str(i)+'[p_dataframe'+str(i)+'$'+columnname+'>='+str(p_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]',file=scriptfile)
            r('p_highlight'+str(i)+'<- p_dataframe'+str(i)+'[p_dataframe'+str(i)+'$'+columnname+'>='+str(p_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]')
            try :
                r('p_highlight'+str(i)+'$geneName<-gsub(pattern="miRNA;",replacement="",x=p_highlight0$geneName)')
                r('p_highlight'+str(i)+'$geneName<-gsub(pattern=";",replacement="",x=p_highlight0$geneName)')
                r('p_highlight'+str(i)+'$geneName<-gsub(pattern="miRNA",replacement="",x=p_highlight0$geneName)')
            except:
                print("no geneName")
            print('p_highlithcolors'+str(i)+'<- rep("red",nrow(p_highlight'+str(i)+'))',file=scriptfile)
            r('p_highlithcolors'+str(i)+'<- rep("red",nrow(p_highlight'+str(i)+'))')
        for i in range(0,len(negtive_winfiles)):
            threshold_title_list=re.split(r"_",negtive_winfiles[i][1].strip())
            n_threshold[i]=float(threshold_title_list[0])
            negtive_filenames[i],negtive_filenameWithPaths[i]=self.prepareMhtFile(negtive_winfiles[i][0], "Hp", "negtive", fillvalue)
            r('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)')
            r('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
            print('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T,stringsAsFactors=FALSE)',file=scriptfile)
            print('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
            r('n_highlight'+str(i)+'<-n_dataframe'+str(i)+'[n_dataframe'+str(i)+'$'+columnname+'<='+str(n_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]')
            try:
                r('n_highlight'+str(i)+'$geneName<-gsub(pattern="miRNA;",replacement="",x=n_highlight0$geneName)')
                r('n_highlight'+str(i)+'$geneName<-gsub(pattern=";",replacement="",x=n_highlight0$geneName)')
                r('n_highlight'+str(i)+'$geneName<-gsub(pattern="miRNA",replacement="",x=n_highlight0$geneName)')
            except:
                print("no geneName")
            print('n_highlight'+str(i)+'<-n_dataframe'+str(i)+'[n_dataframe'+str(i)+'$'+columnname+'<='+str(n_threshold[i])+',][,c("chrNo","winNo","'+columnname+'","geneName")]',file=scriptfile)
            print('n_highlithcolors'+str(i)+'<- rep("red",nrow(n_highlight'+str(i)+'))',file=scriptfile)
            r('n_highlithcolors'+str(i)+'<- rep("red",nrow(n_highlight'+str(i)+'))')
            
        r('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta"),300)')
        r('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(0.8, 4, 0.8, 2))')
        print('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta"),300)',file=scriptfile)
        print('par(las=1, cex.axis=1.5, cex=0.8,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)) +',1),mar=c(0.8, 4, 0.8, 2))',file=scriptfile)
        if len(positive_winfiles)==1:
            hopscex='1'
        else:
            hopscex=str(0.6*(len(positive_winfiles)+len(negtive_winfiles)))
        n=[]
        intersectionSet={}
        intersectionlistMap={}
        for i in range(0,len(positive_winfiles)):
            intersectionlistMap[str(i)+"p"]=[]
            r('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)')
            
            r('hops<-hmht.control(data=p_highlight'+str(i)+',cex='+hopscex+',yoffset=0.2)')
            a=list(r('hops$data$geneName[!is.na(hops$data$geneName)&hops$data$geneName!="top1"]'))
            for e in a:
                nn=re.split(r";",e)
                for ee in nn:
                    n.append(ee)
                    intersectionlistMap[str(i)+"p"].append(ee)
            r('mhtplot(p_data'+str(i)+',ops,hops,pch=19,ylab="z' + "Fst" + '",xlab="")')
            threshold_title_list=re.split(r"_",positive_winfiles[i][1].strip())
            r("title(main='" + " ".join(threshold_title_list[1:]) + "',cex.main=0.8)")
            r('axis(2)')
            r('abline(h='+str(threshold_title_list[0])+')')
            print('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)',file=scriptfile)
            print('hops<-hmht.control(data=p_highlight'+str(i)+',cex='+str(0.6*(len(positive_winfiles)+len(negtive_winfiles)))+',yoffset=0.2)',file=scriptfile)
            print('mhtplot(p_data'+str(i)+',ops,hops,pch=19,ylab="z' + "Fst" + '",xlab="")',file=scriptfile)
            print("title(main='" + positive_winfiles[i][2] + "',cex.main=0.8)",file=scriptfile)
            print('axis(2)',file=scriptfile)
            print('abline(h=5)',file=scriptfile)
        for i in range(0,len(negtive_winfiles)):
            intersectionlistMap[str(i)+"n"]=[]
            r('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)')
            r('hops<-hmht.control(data=n_highlight'+str(i)+',cex='+hopscex+',yoffset=-0.2)')
            a=list(r('hops$data$geneName[!is.na(hops$data$geneName)&hops$data$geneName!="top1"]'))
            for e in a:
                nn=re.split(r";",e)
                for ee in nn:
                    n.append(ee)
                    intersectionlistMap[str(i)+"n"].append(ee)
            print('ops<-mht.control(logscale=FALSE,colors=colors,cex=0.6)',file=scriptfile)
            print('hops<-hmht.control(data=n_highlight'+str(i)+',cex=0.6,yoffset=-0.2)',file=scriptfile)
            print('mhtplot(n_data'+str(i)+',ops,hops,pch=19,ylab="zHp",xlab="")',file=scriptfile)            
            r('mhtplot(n_data'+str(i)+',ops,hops,pch=19,ylab="zHp",xlab="")')
            threshold_title_list=re.split(r"_",negtive_winfiles[i][1].strip())
            r("title(main='" + " ".join(threshold_title_list[1:]) + "',cex.main=0.8)")
            r('axis(2)')
            r('abline(h='+str(threshold_title_list[0])+')')

            print("title(main='" + negtive_winfiles[i][2] + "',cex.main=0.8)",file=scriptfile)
            print('axis(2)',file=scriptfile)
            print('abline(h='+str(threshold_title_list[0])+')',file=scriptfile)
        r('axis(1)')
        r('dev.off()')
        print(r('Cairo.capabilities()'))
        scriptfile.close()
        intersectionSet=intersectionlistMap[list(intersectionlistMap.keys())[0]]
        for k in intersectionlistMap.keys():
            intersectionSet=set(intersectionSet) & set(intersectionlistMap[k])
        return list(set(n)),list(intersectionSet)

    def makeMhtplots_compareInOnePicture(self, outputnamewithpath,positive_winfiles,negtive_winfiles,outfileNameWIN_Alist,fillvalue=0,columnname="zvalue",splitintoparts=1):
        scriptfile=open("stripts.R",'w')
        print(outputnamewithpath,file=scriptfile)
        print(positive_winfiles,negtive_winfiles)
        randomstr=Util.random_str()
        NoOfcurchrom={}
        for winfile,threshold,outbedfilename in positive_winfiles+negtive_winfiles+outfileNameWIN_Alist:
            NoOfcurchrom[winfile]=0
            os.system("awk 'NR>1{print $1}' "+winfile+"|uniq >"+winfile+"_chrom")
        for winfile_idx in range(1,len(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)):
            print("please input winfile withpath")
            a=os.system("comm -12 "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx-1][0]+"_chrom "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx][0]+"_chrom|uniq > temp_chrom"+randomstr)
#             print("comm -12 "+(positive_winfiles+negtive_winfiles)[winfile_idx-1]+"_chrom "+(positive_winfiles+negtive_winfiles)[winfile_idx]+"_chrom|uniq > temp_chrom"+randomstr)
            if a!=0:
                print("comm -12 "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx-1][0]+"_chrom "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx][0]+"_chrom|uniq > temp_chrom"+randomstr)
                exit(-1)   
            print(winfile_idx,"mv temp_chrom"+randomstr+" "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx][0]+"_chrom ")
            os.system("mv temp_chrom"+randomstr+" "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx][0]+"_chrom ")
  
        for winfile_idx in range(1,len(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)):
            os.system("rm "+(positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[winfile_idx-1][0]+"_chrom ")
        print((positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[-1][0]+"_chrom ")
        t=open((positive_winfiles+negtive_winfiles+outfileNameWIN_Alist)[-1][0]+"_chrom","r")
        chromlist=t.readlines();t.close();#os.system("rm "+(positive_winfiles+negtive_winfiles)[-1]+"_chrom ")
        d=len(chromlist)/splitintoparts
        part_i=0
        chromfile=open("chromfile.txt","w")
        for part_i in range(0,splitintoparts-1):
            if re.search(r"^.*/", outputnamewithpath)!=None:
                dir = re.search(r"^.*/", outputnamewithpath).group(0)
            else:
                dir=self.olddir
            outname=re.search(r"[^/]*$", outputnamewithpath).group(0)
            r = robjects.r
            if positive_winfiles!=[]:
                positive_filenames=[""]*len(positive_winfiles);positive_filenameWithPaths=[""]*len(positive_winfiles)
            if negtive_winfiles!=[]:
                negtive_filenames=[""]*len(negtive_winfiles);negtive_filenameWithPaths=[""]*len(negtive_winfiles)
            if outfileNameWIN_Alist!=[]:
                outfileNameWIN_Alistfilenames=[""]*len(outfileNameWIN_Alist);outfileNameWIN_AlistfilenamesWithPaths=[""]*len(outfileNameWIN_Alist)
            print(dir)
            print("setwd('" + dir + "')")
            r("setwd('" + dir + "')")
            r('.libPaths("/home/liurui/software/Rpackages")')
            r("library(gap)")
            print("library(gap)",file=scriptfile)
            r("library(Cairo)")
            print("library(Cairo)",file=scriptfile)
    #         r('Cairo("'+outname+'.pdf",width='+str(((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)*2)+',type="png",height='+str((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)+')')
            templ=positive_winfiles+negtive_winfiles+outfileNameWIN_Alist
#             r('CairoPS("'+outname+chromlist[NoOfcurchrom[templ[0][0]]].strip()+"_"+chromlist[math.ceil(NoOfcurchrom[templ[0][0]]+d)].strip()+"part_"+str(part_i)+'.ps",width=1800,height=900)')
            print('Cairo("'+outname+chromlist[NoOfcurchrom[templ[0][0]]].strip()+"_"+chromlist[math.ceil(NoOfcurchrom[templ[0][0]]+d)].strip()+"part_"+str(part_i)+'.png",type="png",pointsize=12,res=750,width=900,height=550,bg="white",units="mm")',file=scriptfile)
            r('Cairo("'+outname+chromlist[NoOfcurchrom[templ[0][0]]].strip()+"_"+chromlist[math.ceil(NoOfcurchrom[templ[0][0]]+d)].strip()+"part_"+str(part_i)+'.png",type="png",pointsize=12,res=300,width=6000,height=3000)')  
            if outfileNameWIN_Alist!=[]:
                for i in range(0,len(outfileNameWIN_Alist)):
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[outfileNameWIN_Alist[i][0]]].strip()+"""/){start="true"}if($0~/"""+chromlist[math.ceil(NoOfcurchrom[outfileNameWIN_Alist[i][0]]+d)].strip()+"""/){end="true"}if(start=="true" && end!="true"){print $0}}' """+outfileNameWIN_Alist[i][0]+"""> """+outfileNameWIN_Alist[i][0]+"part_"+str(part_i))
                    if i==0:
                        print(NoOfcurchrom[outfileNameWIN_Alist[i][0]],file=chromfile)
                        print(math.ceil(NoOfcurchrom[outfileNameWIN_Alist[i][0]]+d)+1,file=chromfile)
                    NoOfcurchrom[outfileNameWIN_Alist[i][0]]=math.ceil(NoOfcurchrom[outfileNameWIN_Alist[i][0]]+d)+1
                    print(chromlist[NoOfcurchrom[outfileNameWIN_Alist[i][0]]].strip())
                    outfileNameWIN_Alistfilenames[i],outfileNameWIN_AlistfilenamesWithPaths[i]=self.prepareMhtFile(outfileNameWIN_Alist[i][0]+"part_"+str(part_i), "TajimasD", "allvalue", fillvalue)
                    r('a_dataframe'+str(i)+'=read.delim("' + outfileNameWIN_AlistfilenamesWithPaths[i] + '",header=T)')
                    print('a_dataframe'+str(i)+'=read.delim("' + outfileNameWIN_AlistfilenamesWithPaths[i] + '",header=T)',file=scriptfile)
                    r('a_data'+str(i)+' <- with(a_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))')
                    print('a_data'+str(i)+' <- with(a_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))',file=scriptfile)
                    os.system("cp "+outfileNameWIN_AlistfilenamesWithPaths[i]+" "+dir+re.search(r"[^/]*$", outfileNameWIN_AlistfilenamesWithPaths[i]).group(0))
#                     os.system("rm "+outfileNameWIN_AlistfilenamesWithPaths[i]+" "+outfileNameWIN_Alist[i][0]+"part_"+str(part_i))
            if positive_winfiles!=[]:
                for i in range(0,len(positive_winfiles)):
                    print("please input winfile withpath")
                    print("NoOfcurchrom[positive_winfiles[i]]",NoOfcurchrom[positive_winfiles[i][0]],"chromlist[NoOfcurchrom[positive_winfiles[i]]]",chromlist[NoOfcurchrom[positive_winfiles[i][0]]],"d",d,"NoOfcurchrom[positive_winfiles[i]]",NoOfcurchrom[positive_winfiles[i][0]],"chromlist",len(chromlist))
                    print("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[positive_winfiles[i][0]]].strip()+"""/){start="true"}if($0~/"""+chromlist[math.ceil(NoOfcurchrom[positive_winfiles[i][0]]+d)].strip()+"""/){end="true"}if(start=="true" && end!="true"){print $0}}' """+positive_winfiles[i][0]+"""> """+positive_winfiles[i][0]+"part_"+str(part_i))
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[positive_winfiles[i][0]]].strip()+"""/){start="true"}if($0~/"""+chromlist[math.ceil(NoOfcurchrom[positive_winfiles[i][0]]+d)].strip()+"""/){end="true"}if(start=="true" && end!="true"){print $0}}' """+positive_winfiles[i][0]+"""> """+positive_winfiles[i][0]+"part_"+str(part_i))
    #                 os.system("""awk '{if($0~/"""+chromlist[math.ceil(NoOfcurchrom[positive_winfiles[i]]+d*len(chromlist))].strip()+"""/){skip="true"}if(skip!="true"){print $0}}' >"""+positive_winfiles[i]+"part_"+str(part_i))
                    if i==0:
                        print(NoOfcurchrom[positive_winfiles[i][0]],file=chromfile)
                        print(math.ceil(NoOfcurchrom[positive_winfiles[i][0]]+d)+1,file=chromfile)
                    NoOfcurchrom[positive_winfiles[i][0]]=math.ceil(NoOfcurchrom[positive_winfiles[i][0]]+d)+1
                    print(chromlist[NoOfcurchrom[positive_winfiles[i][0]]].strip())
                    positive_filenames[i],positive_filenameWithPaths[i]=self.prepareMhtFile(positive_winfiles[i][0]+"part_"+str(part_i), "Fst", "positive", fillvalue)
                    print('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
                    r('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)')
                    print(positive_filenameWithPaths[i])
                    print('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))',file=scriptfile)
                    r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))')
                    os.system("cp "+positive_filenameWithPaths[i]+" "+dir+re.search(r"[^/]*$", positive_filenameWithPaths[i]).group(0))
#                     os.system("rm "+positive_filenameWithPaths[i]+" "+positive_winfiles[i][0]+"part_"+str(part_i))
    #             r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,zvalue))')
            if negtive_winfiles!=[]:
                for i in range(0,len(negtive_winfiles)):
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[negtive_winfiles[i][0]]].strip()+"""/){start="true"}if($0~/"""+chromlist[math.ceil(NoOfcurchrom[negtive_winfiles[i][0]]+d)].strip()+"""/){end="true"}if(start=="true" && end!="true"){print $0}}' """+negtive_winfiles[i][0]+""" > """+negtive_winfiles[i][0]+"part_"+str(part_i))
    #                 os.system("""awk '{if($0~/"""+chromlist[math.ceil(NoOfcurchrom[negtive_winfiles[i]]+d*len(chromlist))].strip()+"""/){skip="true"}if(skip!="true"){print $0}}' >"""+negtive_winfiles[i]+"part_"+str(part_i))
                    NoOfcurchrom[negtive_winfiles[i][0]]=math.ceil(NoOfcurchrom[negtive_winfiles[i][0]]+d)+1
                    negtive_filenames[i],negtive_filenameWithPaths[i]=self.prepareMhtFile(negtive_winfiles[i][0]+"part_"+str(part_i), "Hp", "negtive", fillvalue)
                    print('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
                    r('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)')
                    print('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
                    r('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
                    os.system("cp "+negtive_filenameWithPaths[i]+" "+dir+re.search(r"[^/]*$", negtive_filenameWithPaths[i]).group(0))
#                     os.system("rm "+negtive_filenameWithPaths[i]+" "+negtive_winfiles[i][0]+"part_"+str(part_i))
            print('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta"),300)',file=scriptfile)
            r('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta"),300)')
            print('par(las=1, cex.axis=1.5, cex=2,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)+len(outfileNameWIN_Alist)) +',1),mar=c(2, 4, 1.5, 2),mgp=c(0,1,-6))',file=scriptfile)
            r('par(las=1, cex.axis=1.5, cex=2,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)+len(outfileNameWIN_Alist)) +',1),mar=c(2, 4, 1.5, 2),mgp=c(0,1,-6))')#rep(0,4)
            if outfileNameWIN_Alist!=[]:
                for i in range(0,len(outfileNameWIN_Alist)):
                    if outfileNameWIN_Alist[i][1].lower().strip()!="none":
                        ylimstr=outfileNameWIN_Alist[i][1].lower().strip().replace("_",",")
                        print(ylimstr)
                        print('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,xlab="",'+"ylim=c("+ylimstr+'))',file=scriptfile)
                        r('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab="",'+"ylim=c("+ylimstr+'))')
                    else:
                        print('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="Tajimas D",xlab="")',file=scriptfile)
                        r('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab="")')
                    if outfileNameWIN_Alist[i][2].lower().strip()!="none":
                        print(outfileNameWIN_Alist[i][2].replace("_"," "))
                        print('title(main=' + outfileNameWIN_Alist[i][2].replace("_"," ") + ',cex.main=2)',file=scriptfile)
                        r('title(main=' + outfileNameWIN_Alist[i][2].replace("_"," ") + ',cex.main=2)')
                    else:
                        print("title(main='" + outfileNameWIN_Alistfilenames[i] + "',cex.main=2)",file=scriptfile)
                        r("title(main='" + outfileNameWIN_Alistfilenames[i] + "',cex.main=2)")
                    print('axis(2,pos=1)')
                    r('axis(2,pos=1)')
            if positive_winfiles!=[]:
                for i in range(0,len(positive_winfiles)):
                    if positive_winfiles[i][1].lower().strip()!="none":
                        ylimstr=positive_winfiles[i][1].lower().strip().replace("_",",")
                        print('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",'  + 'xlab=""'+",ylim=c("+ylimstr+'))',file=scriptfile)
                        r('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",'  + 'xlab=""'+",ylim=c("+ylimstr+'))')
                    else:
                        print('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="z' + "Fst" + '",xlab="")',file=scriptfile)
                        r('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="' + "" + '",xlab="")')
                    if positive_winfiles[i][2].lower().strip()!="none":
                        if positive_winfiles[i][2].find("expression")!=-1:
                            print("title(main=" + positive_winfiles[i][2].replace("_"," ") + ",cex.main=2)",file=scriptfile)
                            print(positive_winfiles[i][2].replace("_"," "))
                            r("title(main=" + positive_winfiles[i][2].replace("_"," ") + ",cex.main=2)")
                        else:
                            print("title(main='" + positive_winfiles[i][2].replace("_"," ") + "',cex.main=2)",file=scriptfile)
                            r("title(main='" + positive_winfiles[i][2].replace("_"," ") + "',cex.main=2)")
                    else:
                        print("title(main='" + positive_filenames[i] + "',cex.main=2)",file=scriptfile)
                        r("title(main='" + positive_filenames[i] + "',cex.main=2)")
                    print('axis(2,pos=1)',file=scriptfile)
                    r('axis(2,pos=1)')
            if negtive_winfiles!=[]:
                for i in range(0,len(negtive_winfiles)):
                    if negtive_winfiles[i][1].lower().strip()!="none":
                        ylimstr=negtive_winfiles[i][1].lower().strip().replace("_",",")
                        print('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab=""'+",ylim=c("+ylimstr+'))',file=scriptfile)
                        r('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab=""'+",ylim=c("+ylimstr+'))')
                    else:
                        print('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="zHp",xlab="")',file=scriptfile)
                        r('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab="")')
                    if negtive_winfiles[i][2].lower().strip()!="none":
                        print("title(main='" + negtive_winfiles[i][2].replace("_"," ") + "',cex.main=2)",file=scriptfile)
                        r("title(main='" + negtive_winfiles[i][2].replace("_"," ") + "',cex.main=2)")
                    else:
                        print("title(main='" + negtive_filenames[i] + "',cex.main=2)",file=scriptfile)
                        r("title(main='" + negtive_filenames[i] + "',cex.main=2)")
                    print('axis(2,pos=1)',file=scriptfile)
                    r('axis(2,pos=1)')
#             r('axis(1)')
            print('dev.off()',file=scriptfile)
            r('dev.off()')
        else:
            chromfile.close()
            part_i+=1
            if re.search(r"^.*/", outputnamewithpath)!=None:
                dir = re.search(r"^.*/", outputnamewithpath).group(0)
            else:
                dir=self.olddir
            outname=re.search(r"[^/]*$", outputnamewithpath).group(0)
            r = robjects.r
            positive_filenames=[""]*len(positive_winfiles);positive_filenameWithPaths=[""]*len(positive_winfiles)
            negtive_filenames=[""]*len(negtive_winfiles);negtive_filenameWithPaths=[""]*len(negtive_winfiles)
            if outfileNameWIN_Alist!=[]:
                outfileNameWIN_Alistfilenames=[""]*len(outfileNameWIN_Alist);outfileNameWIN_AlistfilenamesWithPaths=[""]*len(outfileNameWIN_Alist)
            r("setwd('" + dir + "')")
            r('.libPaths("/home/liurui/software/Rpackages")')
            r("library(gap)")
            r("library(Cairo)")
            templ=positive_winfiles+negtive_winfiles+outfileNameWIN_Alist
    #         r('CairoPNG("'+outname+'.png",width='+str(((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)*2)+',height='+str((len(positive_winfiles)+len(negtive_winfiles))*221.5+35)+')')
#             r('CairoPS("'+outname+chromlist[NoOfcurchrom[templ[0][0]]].strip()+"_"+chromlist[-1].strip()+"part_"+str(part_i)+'.ps")')
            r('Cairo("'+outname+chromlist[NoOfcurchrom[templ[0][0]]].strip()+"_"+chromlist[-1].strip()+"part_"+str(part_i)+'.png",type="png",width=1800,height=900)')
            if outfileNameWIN_Alist!=[]:
                for i in range(0,len(outfileNameWIN_Alist)):
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[outfileNameWIN_Alist[i][0]]].strip()+"""/){start="true"}if(start=="true"){print $0}}' """+ outfileNameWIN_Alist[i][0]+""">"""+outfileNameWIN_Alist[i][0]+"part_"+str(part_i))
                    outfileNameWIN_Alistfilenames[i],outfileNameWIN_AlistfilenamesWithPaths[i]=self.prepareMhtFile(outfileNameWIN_Alist[i][0]+"part_"+str(part_i), "TajimasD", "allvalue", fillvalue)
#                     print('a_dataframe'+str(i)+'=read.delim("' + outfileNameWIN_AlistfilenamesWithPaths[i] + '",header=T)',file=scriptfile)
                    r('a_dataframe'+str(i)+'=read.delim("' + outfileNameWIN_AlistfilenamesWithPaths[i] + '",header=T)')
                    r('a_data'+str(i)+' <- with(a_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))')
#                     print('a_data'+str(i)+' <- a_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
                    os.system("cp "+outfileNameWIN_AlistfilenamesWithPaths[i]+" "+dir+re.search(r"[^/]*$", outfileNameWIN_AlistfilenamesWithPaths[i]).group(0))
#                     os.system("rm "+outfileNameWIN_AlistfilenamesWithPaths[i]+" "+outfileNameWIN_Alist[i][0]+"part_"+str(part_i))
            if positive_winfiles!=[]:
                for i in range(0,len(positive_winfiles)):
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[positive_winfiles[i][0]]].strip()+"""/){start="true"}if(start=="true"){print $0}}' """+ positive_winfiles[i][0]+""">"""+positive_winfiles[i][0]+"part_"+str(part_i))
                    positive_filenames[i],positive_filenameWithPaths[i]=self.prepareMhtFile(positive_winfiles[i][0]+"part_"+str(part_i), "Fst", "positive", fillvalue)
#                     print('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
                    r('p_dataframe'+str(i)+'=read.delim("' + positive_filenameWithPaths[i] + '",header=T)')
                    r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,'+columnname+'))')
#                     print('p_data'+str(i)+' <- p_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
                    os.system("cp "+positive_filenameWithPaths[i]+" "+dir+re.search(r"[^/]*$", positive_filenameWithPaths[i]).group(0))
#                     os.system("rm "+positive_filenameWithPaths[i]+" "+positive_winfiles[i][0]+"part_"+str(part_i))
        #             r('p_data'+str(i)+' <- with(p_dataframe'+str(i)+',cbind(chrNo,winNo,zvalue))')
            if negtive_winfiles!=[]:
                for i in range(0,len(negtive_winfiles)):
                    os.system("""awk '{if(NR==1){print$0}if($0~/"""+chromlist[NoOfcurchrom[negtive_winfiles[i][0]]].strip()+"""/){start="true"}if(start=="true"){print $0}}' """+ negtive_winfiles[i][0]+""">"""+negtive_winfiles[i][0]+"part_"+str(part_i))
                    negtive_filenames[i],negtive_filenameWithPaths[i]=self.prepareMhtFile(negtive_winfiles[i][0]+"part_"+str(part_i), "Hp", "negtive", fillvalue)
                    r('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)')
                    r('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]')
                    os.system("cp "+negtive_filenameWithPaths[i] +" "+dir+re.search(r"[^/]*$", negtive_filenameWithPaths[i]).group(0))
#                     os.system("rm "+negtive_filenameWithPaths[i]+" "+negtive_winfiles[i][0]+"part_"+str(part_i))
#                     print('n_dataframe'+str(i)+'=read.delim("' + negtive_filenameWithPaths[i] + '",header=T)',file=scriptfile)
#                     print('n_data'+str(i)+' <- n_dataframe'+str(i)+'[,c("chrNo","winNo","'+columnname+'")]',file=scriptfile)
            r('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta"),300)')
            r('par(las=1, cex.axis=0.15, cex=0.5,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)+len(outfileNameWIN_Alist)) +',1),mar=rep(0,4))')#,mgp=c(-2,1,-4)
#             print('colors <- rep(c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red"),300)',file=scriptfile)
#             print('par(las=1, cex.axis=2, cex=1,mfrow=c('+str(len(positive_winfiles)+len(negtive_winfiles)+len(outfileNameWIN_Alist)) +',1),mar=c(2, 4, 1.5, 2))',file=scriptfile)
            if outfileNameWIN_Alist!=[]:
                for i in range(0,len(outfileNameWIN_Alist)):
                    if outfileNameWIN_Alist[i][1].lower().strip()!="none":
                        ylimstr=outfileNameWIN_Alist[i][1].lower().strip().replace("_",",")
                        print(ylimstr)
                        r('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab=""'+",ylim=c("+ylimstr+'))')
                    else:
                        r('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab="")')
                    if outfileNameWIN_Alist[i][2].lower().strip()!="none":
                        print(outfileNameWIN_Alist[i][2])
                        r('title(main="' + outfileNameWIN_Alist[i][2].replace("_"," ") + '",cex.main=2)')
                    else:
                        r("title(main='" + outfileNameWIN_Alistfilenames[i] + "',cex.main=2)")
                    r('axis(2,pos=1)')
#                     print('mhtplot(a_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=25,ylab="",ylab="Tajimas D",xlab="")',file=scriptfile)
#                     print("title(main='" + outfileNameWIN_Alistfilenames[i] + "',cex.main=2)",file=scriptfile)
#                     print('axis(2)',file=scriptfile)
            if positive_winfiles!=[]:
                for i in range(0,len(positive_winfiles)):
                    if positive_winfiles[i][1].lower().strip()!="none":
                        ylimstr=positive_winfiles[i][1].lower().strip().replace("_",",")
                        r('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",'  + 'xlab=""'+",ylim=c("+ylimstr+'))')
                    else:
                        r('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="' + "" + '",xlab="")')
                    if positive_winfiles[i][2].lower().strip()!="none":
                        if positive_winfiles[i][2].find("expression")!=-1:
                            print('title(main=' + positive_winfiles[i][2].replace("_"," ") + ',cex.main=2)')
                            r('title(main=' + positive_winfiles[i][2].replace("_"," ") + ',cex.main=2)')
                        else:
                            r("title(main='" + positive_winfiles[i][2].replace("_"," ") + "',cex.main=2)")
                    else:
                        r("title(main='" + positive_filenames[i] + "',cex.main=2)")
                    r('axis(2,pos=1)')
#                     print('mhtplot(p_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=16,' +  '",xlab="")',file=scriptfile)
#                     print("title(main='" + positive_filenames[i] + "',cex.main=2)",file=scriptfile)
#                     print('axis(2)',file=scriptfile)
            if negtive_winfiles!=[]:
                for i in range(0,len(negtive_winfiles)):
                    if negtive_winfiles[i][1].lower().strip()!="none":
                        ylimstr=negtive_winfiles[i][1].lower().strip().replace("_",",")
                        r('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab=""'+",ylim=c("+ylimstr+'))')
                    else:
                        r('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="",xlab="")')
                    if negtive_winfiles[i][2].lower().strip()!="none":
                        if negtive_winfiles[i][2].find("expression")!=-1:
                            r("title(main=" + negtive_winfiles[i][2].replace("_"," ") + ",cex.main=2)")
                        else:
                            r("title(main='" + negtive_winfiles[i][2].replace("_"," ") + "',cex.main=2)")
                    else:
                        r("title(main='" + negtive_filenames[i] + "',cex.main=2)")
                    r('axis(2,pos=1)')
#                     print('mhtplot(n_data'+str(i)+',control=mht.control(logscale=FALSE,colors=colors,cex=0.7),pch=15,ylab="zHp",xlab="")',file=scriptfile)
#                     print("title(main='" + negtive_filenames[i] + "',cex.main=2)",file=scriptfile)
#                     print('axis(2)',file=scriptfile)
#             r('axis(1)')
            r('dev.off()')
            
            
            print(r('Cairo.capabilities()'))
            print(splitintoparts)
#             os.system("rm  "+negtive_filenameWithPaths[i]+" "+negtive_winfiles[i]+"part_"+str(part_i)+"  "+positive_filenameWithPaths[i]+" "+positive_winfiles[i]+"part_"+str(part_i))
        scriptfile.close()
        
