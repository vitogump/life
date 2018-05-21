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
#    print("python CaculateSNPratePerBin.py [vcf1] [vcf2] [vcf3]....-d [dbname] -c [chromtable] -w [winwidth] -s [slidesize] -n [species1] -n [species2] ...")
#    exit(-1)
    
parser = OptionParser()
parser.add_option("-c", "--chromtable", dest="chromtable",# action="callback",type="string",callback=useoptionvalue_previous2,
                  help="write report to FILE")
parser.add_option("-o","--outputpath",dest="outputpath",help="default infile1_infile2")
parser.add_option("-C","--Coveragedbin",dest="coveragebin")
parser.add_option("-I","--howtoIndel",dest="howtoIndel",default="no",help="no just")
parser.add_option("-w","--winwidth",dest="winwidth",help="default infile1_infile2")#
parser.add_option("-n", "--speciesname", action="append",dest="specieses",default=[],help="species name")
parser.add_option("-s","--slidesize",dest="slidesize",help="default infile2_infile1")#
parser.add_option("-m","--minlength",dest="minlength")
parser.add_option("-v","--vcffile",dest="vcffile",action="append", default=[],help="default infile1_infile2")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

howtoIndel=options.howtoIndel.strip()
windowWidth=int(options.winwidth)
slideSize=int(options.slidesize)
outputpath=options.outputpath
chromtable=options.chromtable
minlength=options.minlength
vcffileslist=options.vcffile
sql = "select * from " + chromtable+" where chrlength>="+minlength

class SNPsPerBIN():
    def __init__(self):
        self.SNPsPerBINMap = {}

if __name__ == '__main__':
    
    speicesidxs_inbindepthmap=[]
    if len(vcffileslist[:])==1 and len(options.specieses)!=0 and options.coveragebin!=None:
        bindepth=Util.BinDepth(options.coveragebin)
        for species in options.specieses:
            speicesidxs_inbindepthmap.append(bindepth.speciesname.index("Depth_for_" + species)+2)
        consider_Depth=True
    else:
        consider_Depth=False
    dbtools = dbm.DBTools(Util.ip, Util.username, Util.password, Util.genomeinfodbname)
    print(vcffileslist[:])
    for vcf in vcffileslist[:]:
        vcfname=re.search(r"[^/]*$",vcf).group(0)
        if re.search(r"indvd[^/]+",vcf)!=None:
            snpcounter = Caculators.Caculate_SNPsPerBIN(windowWidth,considerINDEL=howtoIndel,MethodToSeq="indvd")
        elif re.search(r"pool[^/]+",vcf)!=None:
            snpcounter = Caculators.Caculate_SNPsPerBIN(windowWidth,considerINDEL=howtoIndel,MethodToSeq="pool")
        outfile = open(outputpath+vcfname + ".snpperbin"+str(windowWidth)+"_"+str(slideSize), 'w')
        print("chrNo\twinNo\tfirstsnppos\tlastsnppos\tnoofsnps\twinvalue\tzvalue",file=outfile)
        win = Util.Window()
        
        pop = VCFutil.VCF_Data(vcf)  # new a class
        snpbinmap = SNPsPerBIN()
#        pop.getVcfMap(vcf)
        
        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable+" where chrlength>="+minlength)[0][0]
        for i in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(i)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0]
                currentchrLen=int(row[1])
                if currentchrID in pop.VcfIndexMap:
                    vcflist_A_chrom = pop.getVcfListByChrom(currentchrID,MQfilter=None)
                    if currentchrLen==0:currentchrLen=vcflist_A_chrom[-1][0]
                    win.slidWindowOverlap(vcflist_A_chrom, currentchrLen, windowWidth, slideSize, snpcounter)
                    snpbinmap.SNPsPerBINMap[currentchrID]=copy.deepcopy(win.winValueL)
                elif currentchrLen!=0:
                    fillNA=[(0,0,0,0)]
                    for i in range(int((currentchrLen-windowWidth)/slideSize)):
                        fillNA.append((0,0,0,0))
                    snpbinmap.SNPsPerBINMap[currentchrID]=fillNA
        
        winCrossGenome = []
        for chrom in snpbinmap.SNPsPerBINMap.keys():
            for i in range(len(snpbinmap.SNPsPerBINMap[chrom])):
                if snpbinmap.SNPsPerBINMap[chrom][i][3]=='NA':
                    snpbinmap.SNPsPerBINMap[chrom][i]=(0,0,0,0)# replace NA to 0
                winCrossGenome.append(snpbinmap.SNPsPerBINMap[chrom][i][3])
        expectation = numpy.mean(winCrossGenome)
        std0 = numpy.std(winCrossGenome)
        std1 = numpy.std(winCrossGenome, ddof=1)
        del winCrossGenome

        totalChroms = dbtools.operateDB("select","select count(*) from "+chromtable+" where chrlength>="+minlength)[0][0]
        for j in range(0,totalChroms,20):
            currentsql=sql+" order by "+primaryID+" limit "+str(j)+",20"
            result=dbtools.operateDB("select",currentsql)
            for row in result:
                currentchrID=row[0].strip()
                currentchrLen=int(row[1])
                
                if currentchrID in snpbinmap.SNPsPerBINMap:       
        #        for chrom in sorted(hscore.HeterozyMap.keys()):
                    for i in range(len(snpbinmap.SNPsPerBINMap[currentchrID])):
                            if consider_Depth and len(speicesidxs_inbindepthmap)==1:
                                if bindepth.depthbinmap[currentchrID][i][speicesidxs_inbindepthmap[0]]=="filtered":
                                    print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t"+"NA"+"\t" + "NA" + "\t" + 'NA', file=outfile)
                                elif bindepth.depthbinmap[currentchrID][i][speicesidxs_inbindepthmap[0]]=="passed":
                                    try:
                                        log2snpsperbin = numpy.log2(snpbinmap.SNPsPerBINMap[currentchrID][i][2])
                                    except AttributeError:
                                        log2snpsperbin=-99999999999999999
                                    print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][2]) + "\t" +'%.15f'%(snpbinmap.SNPsPerBINMap[currentchrID][i][3])+"\t"+ '%.12f'%(log2snpsperbin), file=outfile)
                                else:
                                    print(currentchrID,str(i),snpbinmap.SNPsPerBINMap[currentchrID],"doesnot exist in snpbinmap")
                            elif consider_Depth and len(speicesidxs_inbindepthmap)>1:# 
                                for idx in speicesidxs_inbindepthmap:# pass if any speicese is "passed"
                                    if bindepth.depthbinmap[currentchrID][i][idx]=="passed":
                                        try:
                                            log2snpsperbin = numpy.log2(snpbinmap.SNPsPerBINMap[currentchrID][i][2])
                                        except AttributeError:
                                            log2snpsperbin=-99999999999999999
                                        print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][2]) + "\t"+'%.15f'%(snpbinmap.SNPsPerBINMap[currentchrID][i][3])+'\t' + '%.12f'%(log2snpsperbin), file=outfile)
                                        break
                                else:
                                    print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) + "\t" +"NA"+"\t"+ "NA" + "\t" + 'NA', file=outfile)
                            elif not consider_Depth:
                                try:
                                    log2snpsperbin = numpy.log2(snpbinmap.SNPsPerBINMap[currentchrID][i][2])
                                except AttributeError:
                                    log2snpsperbin=-99999999999999999
                                print(currentchrID + "\t" + str(i) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][0]) + "\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][1]) +"\t" + str(snpbinmap.SNPsPerBINMap[currentchrID][i][2]) +"\t" + '%.15f'%(snpbinmap.SNPsPerBINMap[currentchrID][i][3]) + "\t" + '%.12f'%(log2snpsperbin), file=outfile)
                else:
                    print("not find this chr:",currentchrID)
        print(vcfname+str(windowWidth)+"_"+str(slideSize), str(expectation), str(std0), str(std1), file=open(outputpath+"caculateSNPratePerBinstaticvalue.txt", 'a'))
        outfile.close()
    dbtools.disconnect()



