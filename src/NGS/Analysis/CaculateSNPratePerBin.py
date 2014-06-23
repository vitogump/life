# -*- coding: UTF-8 -*-
from NGS.BasicUtil import *
from optparse import OptionParser
import NGS.BasicUtil.Util
import re,copy
import numpy
import pickle
import src.NGS.BasicUtil.DBManager as dbm
'''
Created on 2013-7-2

@author: rui
'''

primaryID = "chrID"
#if len(sys.argv) < 4:
#    print("python CaculateSNPratePerBin.py [vcf1] [vcf2] [vcf3]....-d [dbname] -c [chromtable] -w [winwidth] -s [slidesize]")
#    exit(-1)
    
parser = OptionParser()
parser.add_option("-d", "--dbname", dest="dbname",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
# (options, args) = parser.parse_args()
parser.add_option("-C","--Coveragedbin",dest="coveragebin")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-n", "--speciesname", dest="species", help="species name")
parser.add_option("-s","--slidesize",dest="slidesize",help="default infile2_infile1")#
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()


windowWidth=int(options.winwidth)
slideSize=int(options.slidesize)
chromtable=options.chromtable
dbname=options.dbname
sql = "select * from " + chromtable

class SNPsPerBIN():
    def __init__(self):
        self.SNPsPerBINMap = {}

if __name__ == '__main__':
    bindepth=Util.BinDepth(options.coveragebin)
    if len(args[:])==1 and options.species!=None:
        speicesidx_inbindepthmap=bindepth.speciesname.index("Depth_for_" + options.species)+2
        consider_Depth=True
    else:
        consider_Depth=False
    dbtools = dbm.DBTools("10.2.48.96", "root", "1234567", dbname)
    print(args[:])
    for vcf in args[:]:
        
        outfile = open(vcf + ".snpperbin"+str(windowWidth)+"_"+str(slideSize), 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\twinvalue\tzvalue",file=outfile)
        win = Util.Window()
        snpcounter = Caculators.Caculate_SNPsPerBIN()
        pop = VCFutil.VCF_Data(vcf)  # new a class
        snpbinmap = SNPsPerBIN()
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
                    win.slidWindowOverlap(vcflist_A_chrom, currentchrLen, windowWidth, slideSize, snpcounter)
                    snpbinmap.SNPsPerBINMap[currentchrID]=copy.deepcopy(win.winValueL)
                else:
                    fillNA=[(0,0,0)]
                    for i in range(int((currentchrLen-windowWidth)/slideSize)):
                        fillNA.append((0,0,0))
                    snpbinmap.SNPsPerBINMap[currentchrID]=fillNA
        
        winCrossGenome = []
        for chrom in snpbinmap.SNPsPerBINMap.keys():
            for i in range(len(snpbinmap.SNPsPerBINMap[chrom])):
                if snpbinmap.SNPsPerBINMap[chrom][i][2]=='NA':
                    snpbinmap.SNPsPerBINMap[chrom][i]=(0,0,0)# replace NA to 0
                winCrossGenome.append(snpbinmap.SNPsPerBINMap[chrom][i][2])
        expectation = numpy.mean(winCrossGenome)
        std0 = numpy.std(winCrossGenome)
        std1 = numpy.std(winCrossGenome, ddof=1)
        del winCrossGenome

        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable)[0][0]
        for j in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(j)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0].strip()
                currentchrLen=int(row[2])
                
                if currentchrID in snpbinmap.SNPsPerBINMap:       
        #        for chrom in sorted(hscore.HeterozyMap.keys()):
                    for i in range(len(snpbinmap.SNPsPerBINMap[currentchrID])):
                            if consider_Depth:
                                if bindepth.depthbinmap[currentchrID][i][speicesidx_inbindepthmap]=="filtered":
                                    print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" + "NA" + "\t" + 'NA', file=outfile)
                                elif bindepth.depthbinmap[currentchrID][i][speicesidx_inbindepthmap]=="passed":
#                                 snpsperkb = snpbinmap.SNPsPerBINMap[currentchrID][i][2]/(windowWidth/1000)
                                    log2snpsperbin = numpy.log2(snpbinmap.SNPsPerBINMap[currentchrID][i][2])
                                    print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" + '%.18f'%(snpbinmap.SNPsPerBINMap[currentchrID][i][2]) + "\t" + '%.12f'%(log2snpsperbin), file=outfile)
                                else:
                                    print(currentchrID,str(i),snpbinmap.SNPsPerBINMap[currentchrID],"doesnot exist in snpbinmap")
                            else:
                                log2snpsperbin = numpy.log2(snpbinmap.SNPsPerBINMap[currentchrID][i][2])
                                print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" + '%.18f'%(snpbinmap.SNPsPerBINMap[currentchrID][i][2]) + "\t" + '%.12f'%(log2snpsperbin), file=outfile)
                else:
                    print(currentchrID)
        print(vcf, str(expectation), str(std0), str(std1), file=open("caculateSNPratePerBinstaticvalue.txt", 'a'))
        outfile.close()
    dbtools.disconnect()



