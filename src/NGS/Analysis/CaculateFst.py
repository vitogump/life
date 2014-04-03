# -*- coding: UTF-8 -*-
from NGS.BasicUtil import *
from itertools import combinations
import NGS.BasicUtil.Util
import numpy
import pickle
import re
import src.NGS.BasicUtil.DBManager as dbm
import sys
import time
SLEEP_FOR_NEXT_TRY=10

'''
Created on 2013-6-30

@author: rui
'''
if len(sys.argv) < 7:
    print("python CaculateFst.py [vcf1] [vcf2] [vcf3]....[globe_Fst(G)/reletivepaire_Fsts(R)] [winwidth] [slidesize] [chromtable]")
    exit(-1)
windowWidth=int(sys.argv[-3])
slideSize=int(sys.argv[-2])
chromtable = sys.argv[-1]
primaryID = "chrID"

sql = "select * from " + chromtable

class Fst():
    def __init__(self):
        super().__init__()
#        self.doubleVcfMap = {}
        self.FstMapByChrom = {}  # {chr:[(first_snp_pos,last_snp_pos,fst),(),()],chr:[],chr:[]}
        self.distMap = {}
    def alin2PopSnpPos(self, vcfMap1, vcfMap2):
        """input:
        two map fomart like this {chrNo:[(pos,REF,ALT,INFO,FORMAT,sample,...),(pos,REF,ALT,INFO,FORMAT,sample,...),,,,,],chrNo:[],,,,,,}
        output:
        one map like this {chrNo:[(pos,REF,ALT,(INFO,FORMAT,sample,...),(INFO,FORMAT,sample,...)),(,,,(),()),,,,,],chrNo:[],,,}
                                                from pop1                        from pop2
        """
        doubleVcfMap={}
        for currentChrom in vcfMap1.keys():
#             self.FstMapByChrom[currentChrom] = []
            doubleVcfMap[currentChrom] = []

            for SNPrec in vcfMap1[currentChrom]:
                low = 0
                if currentChrom not in vcfMap2:
                    print("alin2PopSnpPos",currentChrom,"didn't find in vcfMap2")
                    break
                high = len(vcfMap2[currentChrom]) - 1
                
                posInPop1 = SNPrec[0]
                RefInPop1 = SNPrec[1]
                AltInPop1 = SNPrec[2]
                if re.search(r"[A-Za-z]+,[A-Za-z]+", AltInPop1) != None:  # multiple allels
                    continue
#                dp4 = re.search(r"DP4=(\d*),(\d*),(\d*),(\d*)", SNPrec[3])
#                 print(dp4.group(0))
                
                while low < high:
                    
                    mid = int((low + high) / 2)
                    if posInPop1 == vcfMap2[currentChrom][mid][0]:
                        if AltInPop1 == vcfMap2[currentChrom][mid][2]:#same alt alle
                            doubleVcfMap[currentChrom].append((posInPop1,RefInPop1,AltInPop1,SNPrec[3:] , vcfMap2[currentChrom][mid][3:]))
                        break
                    elif posInPop1 < vcfMap2[currentChrom][mid][0]:
                        high = mid - 1
                    else:
                        low = mid + 1
                else:
                    pass
#                     self.doubleVcfMap[currentChrom].append(SNPrec+)
        return doubleVcfMap

    def caculateFstAccordingdb(self,dbtools,chromstable,vcfNAME_POP1,vcfNAME_POP2,caculator,winwidth,slideSize):
        pop1 = VCFutil.VCF_Data(vcfNAME_POP1)  # new a class
        pop2 = VCFutil.VCF_Data(vcfNAME_POP2)  # new a class
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromstable)[0][0]
        ########################### caculate Fst across all vcf file and fill in self.FstMapByChrom 
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)

            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[2])
                if currentchrID in pop1.VcfIndexMap:
                    pop1SeqOfAChr={}
                    pop2SeqOfAChr={}
                    pop1SeqOfAChr[currentchrID]=pop1.getVcfListByChrom(vcfNAME_POP1, currentchrID)
                    pop2SeqOfAChr[currentchrID]=pop2.getVcfListByChrom(vcfNAME_POP2, currentchrID)
                    self.caculateFst(pop1SeqOfAChr,pop2SeqOfAChr, fst_caculator,currentchrID,currentchrLen,winwidth,slideSize)
                else:#pop1 don't contation the current chromosome
                    fillNA=[(0,0,'NA')]
                    for i in range(int((currentchrLen-windowWidth)/slideSize)):
                        fillNA.append((0,0,'NA'))
                    self.FstMapByChrom[currentchrID]=fillNA
                                 
    def caculateFst(self, vcfMap1_ref, vcfMap2, caculator,currentchrID,currentchrLen, winwidth, slideSize):
        win = Util.Window()
        tempmap={}
        try:
#            self.doubleVcfMap={}
            doubleVcfMap = self.alin2PopSnpPos(vcfMap1_ref, vcfMap2)#produce self.doubleVcfMap{}
#            for currentChrom in self.doubleVcfMap.keys():
    #             self.FstMapByChrom[currentChrom]=[]
            win.winValueL = []
            print("caculateFst value in "+currentchrID)
            
            win.slidWindowOverlap(doubleVcfMap[currentchrID], currentchrLen,winwidth, slideSize, caculator)
            self.FstMapByChrom[currentchrID] = win.winValueL
        except TypeError:#vcfMap2(pop2) don't contation the current chromosome
            print("caculateFst TypeError")
            fillNA=[(0,0,'NA')]
            for i in range(int((currentchrLen-windowWidth)/slideSize)):
                fillNA.append((0,0,'NA'))
            self.FstMapByChrom[currentchrID]=fillNA           

if __name__ == '__main__':
    dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
    if sys.argv[-4]=='R' or sys.argv[-4]=='r':
        
        allspeices=[]
        tableindextoarrayindex=[]
        treearrayprename=""
        for pathtoname in sys.argv[1:-4]:
            allspeices.append(re.search(r"[^/]*$",pathtoname).group(0).replace('.','_'))
            treearrayprename+=re.search(r"[^/]*$",pathtoname).group(0)[0]
        phyliparrayinfile=open(treearrayprename+"phylip.arrayin"+str(windowWidth)+"_"+str(slideSize),'w')
        print("mysqltablename: "+treearrayprename+"treearray")
        arraytitle=""
        for name in allspeices:
            arraytitle+=(name+"\t")
        print("\t"+arraytitle+"\n")
        for namerow in allspeices:
            print(namerow[0:8]+"\n")        
            
        allkindofpaire = list(combinations(sys.argv[1:-4], 2))
        alldistMap={}
        tempdbtools = dbm.DBTools("localhost", "root", "1234567", "temp")
        TABLES = {}
        TABLES[treearrayprename+"treearray"] = (
            "CREATE TABLE "+treearrayprename+"treearray ("
            " `chrID` varchar(128) NOT NULL ,"
            " `winNo` int(18) NOT NULL,"
            " PRIMARY KEY (`chrID`,`winNo`)"
            ")engine=innodb default charset=utf8"
            )
        tempdbtools.drop_table(treearrayprename+"treearray")
        time.sleep(SLEEP_FOR_NEXT_TRY)
        tempdbtools.create_table(TABLES)        
        for fstpaire in allkindofpaire:

            fstpaire1name = re.search(r"[^/]*$",fstpaire[0]).group(0).replace('.','_')
            fstpaire2name = re.search(r"[^/]*$", fstpaire[1]).group(0).replace('.','_')  # for linux
            tableindextoarrayindex.append((allspeices.index(fstpaire1name),allspeices.index(fstpaire2name)))
            
            outfile = open(fstpaire1name + fstpaire2name + ".fst"+str(windowWidth)+"_"+str(slideSize), 'w')
            
#             win = Util.Window()
            fst_caculator = Caculators.Caculate_Fst()

            fst = Fst() 
            tempdbtools.operateDB("callproc", "mysql_sp_add_column", data=("temp", treearrayprename+"treearray", (fstpaire1name[0:5]+fstpaire2name[0:5]), "text", "default null"))
                
            print("startcaculatefst:\n", fstpaire1name,fstpaire[0],'\n', fstpaire2name,fstpaire[1])
            fst.caculateFstAccordingdb(dbtools, chromtable, fstpaire[0], fstpaire[1], fst_caculator, windowWidth,slideSize)

            winCrossGenome = []
            for chrom in fst.FstMapByChrom.keys():
                for i in range(len(fst.FstMapByChrom[chrom])):
                    if fst.FstMapByChrom[chrom][i][2] != "NA":
                        winCrossGenome.append(fst.FstMapByChrom[chrom][i][2])
            exception = numpy.mean(winCrossGenome)
            std0 = numpy.std(winCrossGenome, ddof=0)
            std1 = numpy.std(winCrossGenome, ddof=1)
            del winCrossGenome
            
            totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
            for i in range(0,totalChroms,20):
                currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
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
                    tempdbtools.operateDB("insert","insert into "+treearrayprename+"treearray(chrID,winNo,"+fstpaire1name[0:5]+fstpaire2name[0:5]+") values(%s,%s,%s) on duplicate key update "+fstpaire1name[0:5]+fstpaire2name[0:5]+" = '"+str(fst.FstMapByChrom[chrom][i][2])+"'",data=(chrom,str(i),str(fst.FstMapByChrom[chrom][i][2])))
            alldistMap[fstpaire1name+fstpaire2name] = (sum / Number,allspeices.index(fstpaire1name))
            outfile.close()
        for n in alldistMap.keys():
            print(n + "\t" + str(alldistMap[n]), file=open("testdist.txt", 'a'))
        tatalwins = tempdbtools.operateDB("select", "select count(*) from "+treearrayprename+"treearray")[0][0]
        for i in range(0, tatalwins, 100):
            wins = tempdbtools.operateDB("select","select * from "+treearrayprename+"treearray order by chrID asc,winNo asc limit "+str(i) +",100")
            for win in wins:
                abandonthisWin=False
                tmparray=[[0 for x in range(len(allspeices))] for y in range(len(allspeices))]
                
#                print("\t"+arraytitle,file=phyliparrayinfile)
                for i in range(len(win[2:])):
                    tmparray[tableindextoarrayindex[i][0]][tableindextoarrayindex[i][1]]=str(win[i+2])
                    tmparray[tableindextoarrayindex[i][1]][tableindextoarrayindex[i][0]]=str(win[i+2])
                    if win[i+2]==None or win[i+2]=="NA" or win[i+2]=='NULL':
                        abandonthisWin=True
                if abandonthisWin:
                    continue
                print("    "+str(len(allspeices)),file=phyliparrayinfile)
                for i in range(len(allspeices)):
                    tmparray[i][i]='0'
                    try:
                        print(allspeices[i][0:8]+"  "+"\t".join(tmparray[i]),file=phyliparrayinfile)
                    except TypeError:
                        print(i,allspeices,tmparray)
                        print('Error:when making phyliparrayinfile')
                        exit(-1)
        tempdbtools.disconnect()
        phyliparrayinfile.close()
    elif sys.argv[-4] == 'G' or sys.argv[-4] == 'g':
        globalFstMapByChrom={}
        fst_caculator = Caculators.Caculate_Fst()

        
#         fst = Fst() 
        specisnum=str(len(sys.argv[1:-4]))
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
                fstlist[-1].caculateFstAccordingdb(dbtools, chromtable, majorpop, othrpop, fst_caculator, windowWidth,slideSize)          
            outfile=open(majorpop+'.gfst'+str(windowWidth)+"_"+str(slideSize)+"_"+specisnum,'w')
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

                totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
                for i in range(0,totalChroms,20):
                    currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
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
            outfile.close()                    
    dbtools.disconnect()
