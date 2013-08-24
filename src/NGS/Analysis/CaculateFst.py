# -*- coding: UTF-8 -*-
from NGS.BasicUtil import *
from itertools import combinations
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm
import re
import numpy
import sys
import pickle


'''
Created on 2013-6-30

@author: rui
'''
if len(sys.argv) < 7:
    print("python CaculateFst.py [vcf1] [vcf2] [vcf3]....[globe_Fst(G)/reletivepaire_Fsts(R)] [winwidth] [slidesize] [fastway]")
    exit(-1)
windowWidth=int(sys.argv[-3])
slideSize=int(sys.argv[-2])
tablename = 'chromosome'
primaryID = "chrID"

sql = "select * from " + tablename

class Fst():
    def __init__(self):
        super().__init__()
        self.doubleVcfMap = {}
        self.FstMapByChrom = {}  # {chr:[(first_snp_pos,last_snp_pos,fst),(),()],chr:[],chr:[]}
        self.distMap = {}
    def alin2PopSnpPos(self, vcfMap1, vcfMap2):
        """
        {chrNo:[(pos,REF,ALT,INFO),(pos,REF,ALT,INFO),,,,,],chrNo:[],,,,,,}
        """
        for currentChrom in vcfMap1.keys():
#             self.FstMapByChrom[currentChrom] = []
            self.doubleVcfMap[currentChrom] = []

            for SNPrec in vcfMap1[currentChrom]:
                low = 0
                if currentChrom not in vcfMap2:
                    break
                high = len(vcfMap2[currentChrom]) - 1
                
                posInPop1 = SNPrec[0]
                RefInPop1 = SNPrec[1]
                AltInPop1 = SNPrec[2]
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
                    continue
                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#                 print(dp4.group(0))
                
                while low < high:
                    
                    mid = int((low + high) / 2)
                    if posInPop1 == vcfMap2[currentChrom][mid][0]:
                        if AltInPop1 == vcfMap2[currentChrom][mid][2]:
                            self.doubleVcfMap[currentChrom].append(SNPrec + vcfMap2[currentChrom][mid])
                        break
                    elif posInPop1 < vcfMap2[currentChrom][mid][0]:
                        high = mid - 1
                    else:
                        low = mid + 1
                else:
                    pass
#                     self.doubleVcfMap[currentChrom].append(SNPrec+)
        return self.doubleVcfMap

    def caculateFstAccordingdb(self,dbtools,chromstable,vcfNAME_POP1,vcfNAME_POP2,caculator,winwidth,slideSize):
        pop1 = VCFutil.VCF_Data(vcfNAME_POP1)  # new a class
        pop2 = VCFutil.VCF_Data(vcfNAME_POP2)  # new a class
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromstable)[0][0]
        ########################### caculate Fst across all vcf file and fill in self.FstMapByChrom 
        for i in range(0,totalChroms,20):
            currentsql=sql+"order by"+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[2])
                if currentchrID in pop1.VcfIndexMap:
                    pop1SeqOfAChr={}
                    pop2SeqOfAChr={}
                    pop1SeqOfAChr[currentchrID]=pop1.getVcfMapByChrom(vcfNAME_POP1, currentchrID)
                    pop2SeqOfAChr[currentchrID]=pop2.getVcfMapByChrom(vcfNAME_POP1, currentchrID)
                    self.caculateFst(pop1SeqOfAChr,pop2SeqOfAChr, fst_caculator,int(sys.argv[-3]),int(sys.argv[-2]))       
    def caculateFst(self, vcfMap1_ref, vcfMap2, caculator, winwidth, slideSize):
        win = Util.Window()
        self.alin2PopSnpPos(vcfMap1_ref, vcfMap2)#produce self.doubleVcfMap{}
        
        for currentChrom in self.doubleVcfMap.keys():
#             self.FstMapByChrom[currentChrom]=[]
            win.winValueL = []
            print("caculateFst value in "+currentChrom)
            win.slidWindowOverlap(self.doubleVcfMap[currentChrom], winwidth, slideSize, caculator)
            self.FstMapByChrom[currentChrom] = win.winValueL


if __name__ == '__main__':
    dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
    if sys.argv[-4]=='R' or sys.argv[-4]=='r':
        allkindofpaire = list(combinations(sys.argv[1:-4], 2))
        alldistMap={}
        for fstpaire in allkindofpaire:
            fstpaire2name = re.search(r"[^/]*$", fstpaire[1]).group(0)  # for linux
            outfile = open(fstpaire[0] + fstpaire2name + ".fst", 'w')
            
#             win = Util.Window()
            fst_caculator = Caculators.Caculate_Fst()

            fst = Fst() 
        
            print("startcaculatefst", fstpaire[0], fstpaire[1])
            fst.caculateFstAccordingdb(dbtools, tablename, fstpaire[1], fstpaire[2], fst_caculator, int(sys.argv[-3]),int(sys.argv[-2]))

            winCrossGenome = []
            for chrom in fst.FstMapByChrom.keys():
                for i in range(len(fst.FstMapByChrom[chrom])):
                    if fst.FstMapByChrom[chrom][i][2] != "NA":
                        winCrossGenome.append(fst.FstMapByChrom[chrom][i][2])
            exception = numpy.mean(winCrossGenome)
            std0 = numpy.std(winCrossGenome, ddof=0)
            std1 = numpy.std(winCrossGenome, ddof=1)
            del winCrossGenome
            
            totalChroms = dbtools.operateDB("select","select count(*) from "+tablename)[0][0]
            for i in range(0,totalChroms,20):
                currentsql=sql+"order by"+primaryID+" limit "+str(i)+",20"
                result=dbtools.operateDB("select",currentsql)
                for row in result:
                    currentchrID=row[0]
                    currentchrLen=int(row[2])
                    if currentchrID in fst.FstMapByChrom:
                        for i in range(len(fst.FstMapByChrom[currentchrID])):
                            if fst.FstMapByChrom[currentchrID][i][2] != "NA":
                                zFst = (fst.FstMapByChrom[currentchrID][i][2] - exception) / std1
                            else:
                                zFst = "NA"
                            print(currentchrID + "\t" + str(i) + "\t" + str(fst.FstMapByChrom[currentchrID][i][0]) + "\t" + str(fst.FstMapByChrom[currentchrID][i][1]) + "\t" + str(fst.FstMapByChrom[currentchrID][i][2]) + "\t" + str(zFst), file=outfile)                        
#            for chrom in sorted(fst.FstMapByChrom.keys()):

            
            sum = 0
            Number = 0
            for chrom in sorted(fst.FstMapByChrom.keys()):
                for i in range(len(fst.FstMapByChrom[chrom])):
                    if fst.FstMapByChrom[chrom][i][2] != 'NA':
                        Number += 1
                        sum += fst.FstMapByChrom[chrom][i][2]
            alldistMap[re.search(r"[^/]*$", fstpaire[0]).group(0) + fstpaire2name] = sum / Number
            outfile.close()
        for n in alldistMap.keys():
            print(n + "\t" + str(alldistMap[n]), file=open("testdist.txt", 'a'))
    elif sys.argv[-4] == 'G' or sys.argv[-4] == 'g':
        globalFstMapByChrom={}
        fst_caculator = Caculators.Caculate_Fst()

        
#         fst = Fst() 
        for majorpop in sys.argv[1:-4]:
#            pop1 = VCFutil.VCF_Data(majorpop)  # new a class
#            pop1.getVcfMap(majorpop)

            fstlist=[]   
            for othrpop in sys.argv[1:-4]:
                if majorpop == othrpop:
                    continue
#                pop2 = VCFutil.VCF_Data(othrpop)  # new a class 
#                pop2.getVcfMap(othrpop)
                print("startcaculatefst", majorpop, othrpop)
                fstlist.append(Fst())
                fstlist[-1].caculateFstAccordingdb(dbtools, tablename, majorpop, othrpop, fst_caculator, int(sys.argv[-3]),int(sys.argv[-2]))
#                fstlist[-1].caculateFst(pop1.VcfMap_AllChrom, pop2.VcfMap_AllChrom, fst_caculator,int(sys.argv[-3]),int(sys.argv[-2]))
            
            outfile=open(majorpop+'.gfst','w')
            if len(fstlist) != 0:
                for chrom in fstlist[0].FstMapByChrom.keys():
                    globalFstMapByChrom[chrom]=[]
                    for winNo in range(0,len(fstlist[0].FstMapByChrom[chrom])):
                        sumFstInAWin=0
                        Number=0
                        for i in range(0,len(fstlist)):
                            try:

                                if fstlist[i].FstMapByChrom[chrom][winNo][2]!= 'NA':
                                    Number+=1
                                    sumFstInAWin+=fstlist[i].FstMapByChrom[chrom][winNo][2]
                            except IndexError:
                                for j in range(0,len(fstlist)):
                                    print(str(j),sys.argv[1+j],chrom,str(winNo),str(len(fstlist[j].FstMapByChrom[chrom])))
                                continue# always in the last position,and the value is caculate any way,so can't mispostion.
                        try:
                            gfst=sumFstInAWin/Number
                        except ZeroDivisionError:
                            gfst="NA"
                        globalFstMapByChrom[chrom].append((fstlist[0].FstMapByChrom[chrom][winNo][0],fstlist[0].FstMapByChrom[chrom][winNo][1],gfst))
#                         print(chrom + "\t" + str(winNo) + "\t" + str(fstlist[0].FstMapByChrom[chrom][winNo][0]) + "\t" + str(fstlist[0].FstMapByChrom[chrom][winNo][1]) + "\t" + str(gfst), file=outfile)


                winCrossGenome = []
                for chrom in globalFstMapByChrom.keys():
                    for i in range(len(globalFstMapByChrom[chrom])):
                        if globalFstMapByChrom[chrom][i][2] != "NA":
                            winCrossGenome.append(globalFstMapByChrom[chrom][i][2])
                exception = numpy.mean(winCrossGenome)
                std0 = numpy.std(winCrossGenome, ddof=0)
                std1 = numpy.std(winCrossGenome, ddof=1)
                del winCrossGenome

                totalChroms = dbtools.operateDB("select","select count(*) from "+tablename)[0][0]
                for i in range(0,totalChroms,20):
                    currentsql=sql+"order by"+primaryID+" limit "+str(i)+",20"
                    result=dbtools.operateDB("select",currentsql)
                    for row in result:
                        currentchrID=row[0]
                        currentchrLen=int(row[2])
                        if currentchrID in globalFstMapByChrom:                                
#                for chrom in sorted(globalFstMapByChrom.keys()):
                            for i in range(len(globalFstMapByChrom[currentchrID])):
                                if globalFstMapByChrom[currentchrID][i][2] != "NA":
                                    zgFst = (globalFstMapByChrom[currentchrID][i][2] - exception) / std1
                                else:
                                    zgFst = "NA"
                                print(currentchrID + "\t" + str(i) + "\t" + str(globalFstMapByChrom[currentchrID][i][0]) + "\t" + str(globalFstMapByChrom[currentchrID][i][1]) + "\t" + str(globalFstMapByChrom[currentchrID][i][2]) + "\t" + str(zgFst), file=outfile)
                                

