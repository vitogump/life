# -*- coding: UTF-8 -*-
import re, numpy, sys, pickle,copy
from NGS.BasicUtil import *
import NGS.BasicUtil.Util
import src.NGS.BasicUtil.DBManager as dbm
from optparse import OptionParser
'''
Created on 2013-7-2

@author: rui
'''

primaryID = "chrID"
parser = OptionParser()
parser.add_option("-d", "--chromdbname", dest="chromdbname",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-o","--outputpath",dest="outputpath",help="default infile1_infile2")
parser.add_option("-v","--vcffile",dest="vcffile",action="append", default=[],help="default infile1_infile2")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-s","--slideSize",dest="slideSize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()




# if len(sys.argv) < 4:
#     print("python CaculateHeterozygosityScore.py [vcf1] [vcf2] [vcf3]....[dbname] [chromtable] [winwidth] [slidesize]")
#     exit(-1)
minlength=options.minlength
windowWidth=int(options.winwidth)
slideSize=int(options.slideSize)
chromtable=options.chromtable
chromdbname=options.chromdbname
outputpath=options.outputpath.strip()
sql = "select * from " + chromtable+" where chrlength>="+minlength
vcffileslist=options.vcffile


class HeterozygosityScore():
    def __init__(self):
        self.HeterozyMap = {}

if __name__ == '__main__':
    dbtools = dbm.DBTools("10.2.48.140", "root", "1234567", chromdbname)
    for vcf in vcffileslist[:]:
        vcfname=re.search(r"[^/]*$",vcf).group(0)
        outfile = open(outputpath+vcfname + ".het"+str(windowWidth)+"_"+str(slideSize), 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnp\twinvalue\tzvalue",file=outfile)
        win = Util.Window()
        if re.search(r"indvd[^/]+",vcf)!=None:
            hp_caculator = Caculators.Caculate_Hp(MethodToSeq="indvd")
        elif re.search(r"pool[^/]+",vcf)!=None:
            hp_caculator = Caculators.Caculate_Hp(MethodToSeq="pool")
        pop = VCFutil.VCF_Data(vcf)  # new a class
        hscore = HeterozygosityScore()
#        pop.getVcfMap(vcf)
        
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable + " where chrlength>="+minlength)[0][0]
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[1])
                if currentchrID in pop.VcfIndexMap:
                    vcflist_A_chrom = pop.getVcfListByChrom(vcf, currentchrID)
                    win.slidWindowOverlap(vcflist_A_chrom, currentchrLen, windowWidth, slideSize, hp_caculator)
                    hscore.HeterozyMap[currentchrID]=copy.deepcopy(win.winValueL)
                else:
                    fillNA=[(0,0,0,'NA')]
                    for i in range(int(currentchrLen/slideSize)):
                        fillNA.append((0,0,0,'NA'))
                    hscore.HeterozyMap[currentchrID]=fillNA
        
        winCrossGenome = []
        for chrom in hscore.HeterozyMap.keys():
            for i in range(len(hscore.HeterozyMap[chrom])):
                if hscore.HeterozyMap[chrom][i][3] != 'NA':
                    winCrossGenome.append(hscore.HeterozyMap[chrom][i][3])
        expectation = numpy.mean(winCrossGenome)
        std0 = numpy.std(winCrossGenome)
        std1 = numpy.std(winCrossGenome, ddof=1)
        del winCrossGenome

        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable+ " where chrlength>="+minlength)[0][0]
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[1])
                if currentchrID in hscore.HeterozyMap:       
        #        for chrom in sorted(hscore.HeterozyMap.keys()):
                    for i in range(len(hscore.HeterozyMap[currentchrID])):
                        if hscore.HeterozyMap[currentchrID][i][3] != "NA":
                            zHp = (hscore.HeterozyMap[currentchrID][i][3] - expectation) / std1
                            print(currentchrID + "\t" + str(i) + "\t" + str(hscore.HeterozyMap[currentchrID][i][0]) + "\t" + str(hscore.HeterozyMap[currentchrID][i][1]) + "\t" +str(hscore.HeterozyMap[currentchrID][i][2])+"\t"+ '%.18f'%(hscore.HeterozyMap[currentchrID][i][3]) + "\t" + '%.12f'%(zHp), file=outfile)
                        else:
                            zHp = "NA"
                            print(currentchrID + "\t" + str(i) + "\t" + str(hscore.HeterozyMap[currentchrID][i][0]) + "\t" + str(hscore.HeterozyMap[currentchrID][i][1]) + "\t"+str(hscore.HeterozyMap[currentchrID][i][2])+"\t" + hscore.HeterozyMap[currentchrID][i][3] + "\t" + zHp, file=outfile)
        print(vcf, str(expectation), str(std0), str(std1), file=open(outputpath+"staticvalue.txt", 'a'))
        outfile.close()
    dbtools.disconnect()



