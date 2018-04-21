# -*- coding: UTF-8 -*-
'''
Created on 2018��4��21��

@author: RuiLiu
'''
import numpy, re, copy,pysam
from optparse import OptionParser

from NGS.BasicUtil import *


parser = OptionParser()
parser.add_option("-v", "--snpfilelist", dest="snpfilelist", action="append",default=[], help="snpfile recode 'chrNo,REF,ALT,P1derFreq,P2derFreq,P3derFreq,P4derFreq,BBBA,ABBA,BABA'")
parser.add_option("-i", "--interval", dest="interval", nargs=3, help="minvalue maxvalue breaks. divid the delta (P1Freq-P2Freq)")
parser.add_option("-a", "--archaicPopConfig", dest="archaic", help="wigeon")
parser.add_option("-o", "--output", dest="output", help="outfileprename")
parser.add_option("-D", "--D_fd_winfile", dest="D_fd_file",default=None, help="D_fd winvalue is D ,zvalue is fd")
(options, args) = parser.parse_args()
minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
dincrease = (maxvalue - minvalue) / breaks
delta_DerAf={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
delta_DerAf.pop(minvalue-dincrease,minvalue);delta_DerAf[(minvalue-dincrease,maxvalue)]={"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
print(delta_DerAf)
WIGEONDEPThreshold=30
if __name__ == '__main__':
    if options.D_fd_file:
        Df=open(options.D_fd_file,"r")
        Df.readline();Dvaluecollector=[]
        for win in Df:
            valuelist=re.split(r"\s+", win.strip())
            if valuelist[5] != "nan":
                Dvaluecollector.append(float(valuelist[5]))
        print(Dvaluecollector)
        meanD=numpy.mean(Dvaluecollector)
        stdD=numpy.std(Dvaluecollector, ddof=1)
        M=len(Dvaluecollector)
        print("len(Dvaluecollector),numy,numpy.std(Dvaluecollector, ddof=1)",M,numpy.std(Dvaluecollector),stdD)
        VarMlist=[]
        for i in range(M):
            Dvaluetemp=copy.deepcopy(Dvaluecollector)
            Dvaluetemp.pop(i)
            iVar=numpy.var(Dvaluetemp)
            VarMlist.append(iVar*M)
        jackknifstd=numpy.std(VarMlist,ddof=1)
        print("D:",meanD,"jackknifstd",jackknifstd)
        print("ZD:",meanD/jackknifstd)
        Df.close()
    #config
    arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE={}
    cf=open(options.archaic,"r")
    for line in cf:
        vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
        if vcffilename_obj!=None:
            vcfname=vcffilename_obj.group(1).strip()
            arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
            arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
        elif line.split():
            arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
    cf.close()
    D_weigonSNPfile=open(options.output+"P123Oweigon.joinSNP","w")
    binfile=open(options.output+".FreqStratifiedBBAA","w")
    #travel
    BBAAcount=BABAcount=ABBAcount=0
    for SnpFile in options.snpfilelist:
        snpfile=open(SnpFile,'r');snpfile.readline();currentchrID=""
        vcflistByChrom={};listOfpopvcfRecsmapByAChr=[]
        for line in snpfile:
            linelist = re.split(r'\s+', line.strip())
            chrom = linelist[0].strip()
            pos = int(linelist[1].strip())
            anc = linelist[2].strip()
            der = linelist[3].strip()
            p1 = float(linelist[4].strip())#mallard population
            p2 = float(linelist[5].strip())#spot-billed population
            for a,b in sorted(delta_DerAf.keys()):
                if (p1-p2)>a and (p1-p2)<=b:
                    if p1>0 or p2>0:
                        delta_DerAf[(a,b)]["BinP1orP2"].append(linelist)
                    if linelist[-3]> linelist[-2] and linelist[-3]>linelist[-1]:
                        delta_DerAf[(a,b)]["BBAA"].append(linelist);BBAAcount+=1
                    elif linelist[-2]>linelist[-3] and linelist[-2]>linelist[-1]:
                        delta_DerAf[(a,b)]["ABBA"].append(linelist);ABBAcount+=1
                    elif linelist[-1]>linelist[-3] and linelist[-1]>linelist[-2]:
                        delta_DerAf[(a,b)]["BABA"].append(linelist);BABAcount+=1
                break
            if chrom!=currentchrID and currentchrID in vcflistByChrom:
                listOfpopvcfRecsmapByAChr.append(vcflistByChrom,{currentchrID:arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].getVcfListByChrom(currentchrID)})
                target_ref_SNPs = Util.alinmultPopSnpPos(listOfpopvcfRecsmapByAChr, "o")
                for cc in target_ref_SNPs.keys():
                    for T in target_ref_SNPs[cc]:
                        if T[3]==None:
                            continue
                        for poprec in T[4:5]:
                            if poprec==None:
                                sum_depth=0
                                for samfile in arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][1:]:
                                    ACGTdep=samfile.count_coverage(currentchrID,T[0]-1,T[0])
                                    for dep in ACGTdep:
                                        sum_depth+=dep[0]
                                if sum_depth>WIGEONDEPThreshold:
                                    AF=0
                                else:
                                    AF="unknow"
                        print(cc,*T[:],AF,file=D_weigonSNPfile)
                currentchrID=chrom
                vcflistByChrom={currentchrID:[]}
            else:
                vcflistByChrom[currentchrID].append((pos,anc,der,linelist[4:]))
    print("BBAAcount,ABBAcount,BABAcount",BBAAcount,ABBAcount,BABAcount)
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(delta_DerAf[(a,b)]):
        print(k,end="\t",file=binfile)
    print("",file=binfile)
    for a,b in sorted(delta_DerAf.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAf[(a,b)]):
            print(len(delta_DerAf[(a,b)][k]),end="\t",file=binfile)
        print("",file=binfile)
    binfile.close();D_weigonSNPfile.close()