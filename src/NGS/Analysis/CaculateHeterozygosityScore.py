# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle
from NGS.BasicUtil import *
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm
'''
Created on 2013-7-2

@author: rui
'''

primaryID = "chrID"



if len(sys.argv) < 4:
    print("python CaculateHeterozygosityScore.py [vcf1] [vcf2] [vcf3]....[dbname] [chromtable] [winwidth] [slidesize]")
    exit(-1)
windowWidth=int(sys.argv[-2])
slideSize=int(sys.argv[-1])
chromtable=sys.argv[-3]
dbname=sys.argv[-4]
sql = "select * from " + chromtable

class HeterozygosityScore():
    def __init__(self):
        self.HeterozyMap = {}

if __name__ == '__main__':
    dbtools = dbm.DBTools("localhost", "root", "1234567", dbname)
    for vcf in sys.argv[1:-4]:
        outfile = open(vcf + ".het"+str(windowWidth)+"_"+str(slideSize), 'w')
        win = Util.Window()
        hp_caculator = Caculators.Caculate_Hp()
        pop = VCFutil.VCF_Data(vcf)  # new a class
        hscore = HeterozygosityScore()
#        pop.getVcfMap(vcf)
        
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[2])
                if currentchrID in pop.VcfIndexMap:
                    vcflist_A_chrom = pop.getVcfListByChrom(vcf, currentchrID)
                    win.slidWindowOverlap(vcflist_A_chrom, currentchrLen, windowWidth, slideSize, hp_caculator)
                    hscore.HeterozyMap[currentchrID]=win.winValueL
                else:
                    fillNA=[(0,0,'NA')]
                    for i in range(int((currentchrLen-windowWidth)/slideSize)):
                        fillNA.append((0,0,'NA'))
                    hscore.HeterozyMap[currentchrID]=fillNA
        
        winCrossGenome = []
        for chrom in hscore.HeterozyMap.keys():
            for i in range(len(hscore.HeterozyMap[chrom])):
                if hscore.HeterozyMap[chrom][i][2] != 'NA':
                    winCrossGenome.append(hscore.HeterozyMap[chrom][i][2])
        expectation = numpy.mean(winCrossGenome)
        std0 = numpy.std(winCrossGenome)
        std1 = numpy.std(winCrossGenome, ddof=1)
        del winCrossGenome

        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[2])
                if currentchrID in hscore.HeterozyMap:       
        #        for chrom in sorted(hscore.HeterozyMap.keys()):
                    for i in range(len(hscore.HeterozyMap[currentchrID])):
                        if hscore.HeterozyMap[currentchrID][i][2] != "NA":
                            zHp = (hscore.HeterozyMap[currentchrID][i][2] - expectation) / std1
                            print(currentchrID + "\t" + str(i) + "\t" + str(hscore.HeterozyMap[currentchrID][i][0]) + "\t" + str(hscore.HeterozyMap[currentchrID][i][1]) + "\t" + '%.18f'%(hscore.HeterozyMap[currentchrID][i][2]) + "\t" + '%.12f'%(zHp), file=outfile)
                        else:
                            zHp = "NA"
                            print(currentchrID + "\t" + str(i) + "\t" + str(hscore.HeterozyMap[currentchrID][i][0]) + "\t" + str(hscore.HeterozyMap[currentchrID][i][1]) + "\t" + hscore.HeterozyMap[currentchrID][i][2] + "\t" + zHp, file=outfile)
        print(vcf, str(expectation), str(std0), str(std1), file=open("staticvalue.txt", 'a'))
        outfile.close()
    dbtools.disconnect()



