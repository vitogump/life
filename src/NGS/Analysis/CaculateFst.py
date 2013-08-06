# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle
from NGS.BasicUtil import *
import NGS.BasicUtil.Util


from itertools import combinations
'''
Created on 2013-6-30

@author: rui
'''
if len(sys.argv) < 7:
    print("python CaculateFst.py [vcf1] [vcf2] [vcf3]....[globe_Fst(G)/reletivepaire_Fsts(R)] [winwidth] [slidesize] [fastway]")
    exit(-1)

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
    if sys.argv[-4]=='R' or sys.argv[-4]=='r':
        allkindofpaire = list(combinations(sys.argv[1:-4], 2))
        alldistMap={}
        for fstpaire in allkindofpaire:
            fstpaire2name = re.search(r"[^/]*$", fstpaire[1]).group(0)  # for linux
            outfile = open(fstpaire[0] + fstpaire2name + ".fst", 'w')
            
#             win = Util.Window()
            fst_caculator = Caculators.Caculate_Fst()
    
            pop1 = VCFutil.VCF_Data()  # new a class
            pop2 = VCFutil.VCF_Data()  # new a class    
    
            fst = Fst() 
        
            if sys.argv[-1] == "slowway":
                try:
                    vcf_1_idx = pickle.load(open(fstpaire[0] + ".myindex", 'rb'))
                    vcf_2_idx = pickle.load(open(fstpaire[1] + ".myindex", 'rb'))
                    
                except IOError:
                    pop1.indexVCF(fstpaire[0], fstpaire[0] + ".myindex")
                    pop2.indexVCF(fstpaire[1], fstpaire[1] + ".myindex")
                    vcf_1_idx = pickle.load(open(fstpaire[0] + ".myindex", 'rb'))
                    vcf_2_idx = pickle.load(open(fstpaire[1] + ".myindex", 'rb'))
                tmppopmap1 = {}
                tmppopmap2 = {}            
                for chrom in vcf_1_idx.keys():
                    if chrom == "title":
                        continue
                    pop1.getVcfMapByChrom(fstpaire[0], chrom, vcf_1_idx)
                    if pop2.getVcfMapByChrom(fstpaire[1], chrom, vcf_2_idx) == -1:
                        continue
                    tmppopmap1[chrom] = pop1.VcfList_A_Chrom
                    tmppopmap2[chrom] = pop2.VcfList_A_Chrom
                    fst.caculateFst(tmppopmap1, tmppopmap2, fst_caculator,int(sys.argv[-3]),int(sys.argv[-2]))
                    for e in fst.FstMapByChrom[chrom]:
                        print(chrom, e[0], e[1], e[2], sep='\t', file=outfile)
                    del tmppopmap1[chrom]
                    del tmppopmap2[chrom]
            elif sys.argv[-1] == "fastway":
                pop1.getVcfMap(fstpaire[0])
                pop2.getVcfMap(fstpaire[1])
                print("startcaculatefst", fstpaire[0], fstpaire[1])
                fst.caculateFst(pop1.VcfMap_AllChrom, pop2.VcfMap_AllChrom, fst_caculator,int(sys.argv[-3]),int(sys.argv[-2]))
    #             for chrom in fst.FstMapByChrom.keys():
    #                 for e in fst.FstMapByChrom[chrom]:
    #                     print(chrom,e[0],e[1],e[2],sep='\t',file=outfile)
                winCrossGenome = []
                for chrom in fst.FstMapByChrom.keys():
                    for i in range(len(fst.FstMapByChrom[chrom])):
                        if fst.FstMapByChrom[chrom][i][2] != "NA":
                            winCrossGenome.append(fst.FstMapByChrom[chrom][i][2])
                exception = numpy.mean(winCrossGenome)
                std0 = numpy.std(winCrossGenome, ddof=0)
                std1 = numpy.std(winCrossGenome, ddof=1)
                del winCrossGenome
                for chrom in sorted(fst.FstMapByChrom.keys()):
                    for i in range(len(fst.FstMapByChrom[chrom])):
                        if fst.FstMapByChrom[chrom][i][2] != "NA":
                            zFst = (fst.FstMapByChrom[chrom][i][2] - exception) / std1
                        else:
                            zFst = "NA"
                        print(chrom + "\t" + str(i) + "\t" + str(fst.FstMapByChrom[chrom][i][0]) + "\t" + str(fst.FstMapByChrom[chrom][i][1]) + "\t" + str(fst.FstMapByChrom[chrom][i][2]) + "\t" + str(zFst), file=outfile)
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
            pop1 = VCFutil.VCF_Data()  # new a class
            pop1.getVcfMap(majorpop)

            fstlist=[]   
#             outfile=open(majorpop+'.gfst','w')
#             if len(fstlist) != 0:
#                 for chrom in fstlist[0].FstMapByChrom.keys():
#                     for winNo in fstlist[0].FstMapByChrom[chrom]:
#                         sumFstInAWin=0
#                         Number=0
#                         for i in fstlist:
#                             if fstlist[0].FstMapByChrom[chrom][winNo][0] != fstlist[i].FstMapByChrom[chrom][winNo][0] or fstlist[0].FstMapByChrom[chrom][winNo][1] != fstlist[i].FstMapByChrom[chrom][winNo][1]:
#                                 print(majorpop+"de shang yi ge"+chrom+)
#                                 exit(-1)
#                             if fstlist[i].FstMapByChrom[chrom][winNo]!= 'NA':
#                                 Number+=1
#                                 sumFstInAWin+=fstlist[i].FstMapByChrom[chrom][winNo]
#                         gfst=sumFstInAWin/Number
#                         print(chrom + "\t" + str(winNo) + "\t" + str(fstlist[0].FstMapByChrom[chrom][winNo][0]) + "\t" + str(fstlist[0].FstMapByChrom[chrom][winNo][1]) + "\t" + str(gfst), file=outfile)
#                 fstlist=[]
            for othrpop in sys.argv[1:-4]:
                if majorpop == othrpop:
                    continue
                pop2 = VCFutil.VCF_Data()  # new a class 
                pop2.getVcfMap(othrpop)
                print("startcaculatefst", majorpop, othrpop)
                fstlist.append(Fst())
                fstlist[-1].caculateFst(pop1.VcfMap_AllChrom, pop2.VcfMap_AllChrom, fst_caculator,int(sys.argv[-3]),int(sys.argv[-2]))
            
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
                for chrom in sorted(globalFstMapByChrom.keys()):
                    for i in range(len(globalFstMapByChrom[chrom])):
                        if globalFstMapByChrom[chrom][i][2] != "NA":
                            zgFst = (globalFstMapByChrom[chrom][i][2] - exception) / std1
                        else:
                            zgFst = "NA"
                        print(chrom + "\t" + str(i) + "\t" + str(globalFstMapByChrom[chrom][i][0]) + "\t" + str(globalFstMapByChrom[chrom][i][1]) + "\t" + str(globalFstMapByChrom[chrom][i][2]) + "\t" + str(zgFst), file=outfile)
                        

