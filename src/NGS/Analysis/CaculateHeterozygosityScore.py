# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle
from NGS.BasicUtil import *
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm
'''
Created on 2013-7-2

@author: rui
'''
tablename = 'chromosome'
primaryID = "chrID"

sql = "select * from " + tablename

if len(sys.argv) < 4:
    print("python CaculateHeterozygosityScore.py [vcf1] [vcf2] [vcf3]....[winwidth] [slidesize]")
    exit(-1)

class HeterozygosityScore():
    def __init__(self):
        self.HeterozyMap = {}

if __name__ == '__main__':
    for vcf in sys.argv[1:-2]:
        outfile = open(vcf + ".het", 'w')
        win = Util.Window()
        hp_caculator = Caculators.Caculate_Hp()
        pop = VCFutil.VCF_Data()  # new a class
        hscore = HeterozygosityScore()
        pop.getVcfMap(vcf)
        dbtools = dbm.DBTools("localhost", "root", "1234567", "life_pilot")
        for chrom in pop.VcfMap_AllChrom.keys():
            win.slidWindowOverlap(pop.VcfMap_AllChrom[chrom], int(sys.argv[-2]), int(sys.argv[-1]), hp_caculator)
            hscore.HeterozyMap[chrom] = win.winValueL
        winCrossGenome = []
        for chrom in hscore.HeterozyMap.keys():
            for i in range(len(hscore.HeterozyMap[chrom])):
                if hscore.HeterozyMap[chrom][i][2] != 'NA':
                    winCrossGenome.append(hscore.HeterozyMap[chrom][i][2])
        expectation = numpy.mean(winCrossGenome)
        std0 = numpy.std(winCrossGenome)
        std1 = numpy.std(winCrossGenome, ddof=1)
        del winCrossGenome
        for chrom in sorted(hscore.HeterozyMap.keys()):
            for i in range(len(hscore.HeterozyMap[chrom])):
                if hscore.HeterozyMap[chrom][i][2] != "NA":
                    zHp = (hscore.HeterozyMap[chrom][i][2] - expectation) / std1
                else:
                    zHp = "NA"
                print(chrom + "\t" + str(i) + "\t" + str(hscore.HeterozyMap[chrom][i][0]) + "\t" + str(hscore.HeterozyMap[chrom][i][1]) + "\t" + str(hscore.HeterozyMap[chrom][i][2]) + "\t" + str(zHp), file=outfile)
        print(vcf, str(expectation), str(std0), str(std1), file=open("staticvalue.txt", 'a'))
        outfile.close()



