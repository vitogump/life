# -*- coding: UTF-8 -*-
'''
Created on 2018��4��21��

@author: RuiLiu
'''
import numpy, re, copy,pysam,pickle,os,time
from optparse import OptionParser
from multiprocessing.dummy import Pool
from src.NGS.BasicUtil import *


parser = OptionParser()
parser.add_option("-v", "--snpfilelist", dest="snpfilelist", action="append",default=[], help="snpfile recode 'chrNo,REF,ALT,P1derFreq,P2derFreq,P3derFreq,P4derFreq,BBBA,ABBA,BABA'")
parser.add_option("-l", "--filelist", dest="filelistfile",default=None,  help="snpfile recode 'chrNo,REF,ALT,P1derFreq,P2derFreq,P3derFreq,P4derFreq,BBBA,ABBA,BABA'")
parser.add_option("-i", "--interval", dest="interval", nargs=3, help="minvalue maxvalue breaks. divid the delta (P1Freq-P2Freq)")
parser.add_option("-a", "--archaicPopConfig", dest="archaic", default=None,help="wigeon")
parser.add_option("-o", "--output", dest="output", help="outfileprename")
parser.add_option("-D", "--D_fd_winfile", dest="D_fd_file",default=None, help="D_fd winvalue is D ,zvalue is fd")
(options, args) = parser.parse_args()
minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
dincrease = (maxvalue - minvalue) / breaks
delta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    delta_DerAftotal[(minvalue,minvalue + dincrease)]={"BinP1andP2":0,"BinP1orP2":0,"BBAA":0,"BABA":0,"ABBA":0,"type1meanvalue":0,"type2meanvalue":0}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
delta_DerAftotal.pop(minvalue-dincrease,minvalue);delta_DerAftotal[(minvalue-dincrease,maxvalue)]={"BinP1andP2":0,"BinP1orP2":0,"BBAA":0,"BABA":0,"ABBA":0,"type1meanvalue":0,"type2meanvalue":0}
delta_DerAftotalP1P3=copy.deepcopy(delta_DerAftotal)
delta_DerAftotalP2P3=copy.deepcopy(delta_DerAftotal)
# print(delta_DerAf)
WIGEONDEPThreshold=30
if options.filelistfile==None:
    filelistfile=open(options.output+'prefix.filelist','w')
    filelistfile.close()
BBAAcountsum=ABBAcountsum=BABAcountsum=0
def countJoinFile(FileNamrPre):
    FileNamrPre=FileNamrPre.strip()
    delta_DerAf={}
    minvalue =float(options.interval[0])
    maxvalue= float(options.interval[1])
    breaks= int(options.interval[2])
    dincrease= (maxvalue - minvalue) / breaks
    while minvalue<=maxvalue - dincrease :
        print(minvalue,minvalue + dincrease)
        delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[],"type1meanvalue":[],"type2meanvalue":[]}
        minvalue += dincrease
    delta_DerAf.pop(minvalue-dincrease,minvalue);delta_DerAf[(minvalue-dincrease,maxvalue)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[],"type1meanvalue":[],"type2meanvalue":[]}
    print(delta_DerAf)
    delta_DerAfP1P3=copy.deepcopy(delta_DerAf)
    delta_DerAfP2P3=copy.deepcopy(delta_DerAf)
    BBAAcount=BABAcount=ABBAcount=0
    snpfile=open(FileNamrPre+".snp",'r');print(snpfile.readline())
    for line in snpfile:
        linelist = re.split(r'\s+', line.strip())
#         print(linelist)
#         return 
        chrom = linelist[0].strip()
        pos = int(linelist[1].strip())
        anc = linelist[2].strip()
        der = linelist[3].strip()
        p=[0,1,2,3]
        p[1] = float(linelist[4].strip())#mallard population
        p[2] = float(linelist[5].strip())#spot-billed population
        p[3] = float(linelist[6].strip())#domestic population
        i=1;j=2;k=0
        for dDAF in [delta_DerAf,delta_DerAfP2P3,delta_DerAfP1P3]:
            k+=1
            if k==1:
                i=1;j=2#P1P2 FOCUS ON BBAA
            elif k==2:
                i=2;j=3#P2P3 FOCUS ON ABBA
            elif k==3:
                i=1;j=3#P1P3 FOCUS ON BABA
            for a,b in sorted(dDAF.keys()):
    #                 print("\t",(p1-p2),">",a,":",(p1-p2)>a,(p1-p2),"<",b,":",(p1-p2)<=b)
                if (p[i]-p[j])>=a and (p[i]-p[j])<=b:
    #                     print("\t","passed bin threshold",p1,p2,linelist[-3],linelist[-2],linelist[-1])
                    dDAF[(a,b)]["type1meanvalue"].append(float(linelist[-k]))
                    if p[i]>0 and p[j]>0:
                        dDAF[(a,b)]["BinP1andP2"].append(linelist)
                    if p[i]>0 or p[j]>0:
                        dDAF[(a,b)]["BinP1orP2"].append(linelist)
    #                         print("\t","BinP1orP2",p1,p2)
                    judge=0
                    if float(linelist[8])> float(linelist[9]) and float(linelist[8])>float(linelist[10]):
    #                         print("\t","BBAA",linelist[-3])
                        dDAF[(a,b)]["BBAA"].append(linelist[-3]);BBAAcount+=1;judge=1
                    elif float(linelist[9])>float(linelist[8]) and float(linelist[9])>float(linelist[10]):
                        dDAF[(a,b)]["ABBA"].append(linelist[-2]);ABBAcount+=1;judge=2
    #                         print("\t","ABBA",linelist[-2])
                    elif float(linelist[10])>float(linelist[8]) and float(linelist[10])>float(linelist[9]):
                        dDAF[(a,b)]["BABA"].append(linelist[-1]);BABAcount+=1;judge=3
    #                         print("\t","BABA",linelist[-1])
                    if judge==k:
                        dDAF[(a,b)]["type2meanvalue"].append(float(linelist[-k]))
                    break
    pickle.dump(delta_DerAf,open(FileNamrPre+".FreqStratifiedP1P2BBAA", 'wb'))
    pickle.dump(delta_DerAfP1P3,open(FileNamrPre+".FreqStratifiedP1P3BABA", 'wb'))
    pickle.dump(delta_DerAfP2P3,open(FileNamrPre+".FreqStratifiedP2P3ABBA", 'wb'))
    print("BBAAcount,ABBAcount,BABAcount",BBAAcount,ABBAcount,BABAcount)
    global BBAAcountsum
    global ABBAcountsum
    global BABAcountsum
    BBAAcountsum+=BBAAcount;ABBAcountsum+=ABBAcount;BABAcountsum+=BABAcount
def travelCountjoinWg(SnpFile):
    #config
    delta_DerAf={}
    minvalue =float(options.interval[0])
    maxvalue= float(options.interval[1])
    breaks= int(options.interval[2])
    dincrease= (maxvalue - minvalue) / breaks
    while minvalue<=maxvalue - dincrease :
        print(minvalue,minvalue + dincrease)
        delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
        minvalue += dincrease
    delta_DerAf.pop(minvalue-dincrease,minvalue);delta_DerAf[(minvalue-dincrease,maxvalue)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    print(delta_DerAf)
    rstr=Util.random_str()
    arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE={}
    if options.archaic !=None:
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
        D_weigonSNPfile=open(options.output+rstr+"P123Oweigon.joinSNP","w")
    filelistfile=open(options.output+'prefix.filelist','a')
    print(options.output+rstr,file=filelistfile)
    filelistfile.close()
    
    
#         binfile=open(options.output+rstr+".FreqStratifiedBBAA","w")
    
    BBAAcount=BABAcount=ABBAcount=0
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
#             print("snpline",linelist,"delta",(p1-p2))
        for a,b in sorted(delta_DerAf.keys()):
#                 print("\t",(p1-p2),">",a,":",(p1-p2)>a,(p1-p2),"<",b,":",(p1-p2)<=b)
            if (p1-p2)>a and (p1-p2)<=b:
#                     print("\t","passed bin threshold",p1,p2,linelist[-3],linelist[-2],linelist[-1])
                if p1>0 and p2>0:
                    delta_DerAf[(a,b)]["BinP1andP2"].append(linelist)
                if p1>0 or p2>0:
                    delta_DerAf[(a,b)]["BinP1orP2"].append(linelist)
#                         print("\t","BinP1orP2",p1,p2)
                if float(linelist[8])> float(linelist[9]) and float(linelist[8])>float(linelist[10]):
#                         print("\t","BBAA",linelist[-3])
                    delta_DerAf[(a,b)]["BBAA"].append(linelist[-3]);BBAAcount+=1
                elif float(linelist[9])>float(linelist[8]) and float(linelist[9])>float(linelist[10]):
                    delta_DerAf[(a,b)]["ABBA"].append(linelist[-2]);ABBAcount+=1
#                         print("\t","ABBA",linelist[-2])
                elif float(linelist[10])>float(linelist[8]) and float(linelist[10])>float(linelist[9]):
                    delta_DerAf[(a,b)]["BABA"].append(linelist[-1]);BABAcount+=1
#                         print("\t","BABA",linelist[-1])
                break
        if options.archaic !=None and chrom!=currentchrID and currentchrID in vcflistByChrom:
            listOfpopvcfRecsmapByAChr=[vcflistByChrom,{currentchrID:arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname][0].getVcfListByChrom(currentchrID)}]
            target_ref_SNPs = Util.alinmultPopSnpPos(listOfpopvcfRecsmapByAChr, "l") 
            for cc in target_ref_SNPs.keys():
                for T in target_ref_SNPs[cc]:
                    if T[3]==None:
                        print("exit only in wigeon")
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
                        else:
                            AF=float(re.search(r"AF=([\d\.e-]+)[;,]", T[4][0]).group(1))
                    if "," in T[2] and re.split(r",",T[2])[1][0] == re.split(r",",T[2])[0] and re.split(r",",T[2])[1][1] == T[1]:
                        print(cc,*T[:3],*(T[3][0]),(1-AF),sep="\t",file=D_weigonSNPfile)
                    elif "," not in T[2] :
                        print(cc,*T[:3],*(T[3][0]),AF,sep="\t",file=D_weigonSNPfile)
                    else:
                        print("exclude dif SNP",T)
            currentchrID=chrom
            vcflistByChrom={currentchrID:[(pos,anc,der,linelist[4:])]}
        elif options.archaic !=None :
            if currentchrID=="":#first line
                currentchrID=chrom;vcflistByChrom={currentchrID:[]}
            vcflistByChrom[currentchrID].append((pos,anc,der,linelist[4:]))
    pickle.dump(delta_DerAf,open(options.output+rstr+".FreqStratifiedBBAA", 'wb'))
    print("BBAAcount,ABBAcount,BABAcount",BBAAcount,ABBAcount,BABAcount)
    global BBAAcountsum
    global ABBAcountsum
    global BABAcountsum
    BBAAcountsum+=BBAAcount;ABBAcountsum+=ABBAcount;BABAcountsum+=BABAcount
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
#     arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE={}
#     cf=open(options.archaic,"r")
#     for line in cf:
#         vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
#         if vcffilename_obj!=None:
#             vcfname=vcffilename_obj.group(1).strip()
#             arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
#             arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(VCFutil.VCF_Data(vcfname))
#         elif line.split():
#             arcpopvcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.Samfile(line.strip(),'rb'))
#     cf.close()
#     D_weigonSNPfile=open(options.output+"P123Oweigon.joinSNP","w")
#     binfile=open(options.output+".FreqStratifiedBBAA","w")
#     print("chrNo\tpos\tANC\tDER\tP1derFreq\tP2derFreq\tP3derFreq\tP4derFreq\tBBBA\tABBA\tBABA\twigeonAF",file=D_weigonSNPfile)
    #travel

    if options.filelistfile!=None:
        filelistfile=open(options.filelistfile,'r')
        flist=filelistfile.readlines()
        pool=Pool(int(len(flist)))
         
        pool.map(countJoinFile,flist)
        pool.close()
        pool.join()
        time.sleep(60)
    else:
        pool=Pool(int(len(options.snpfilelist)))
         
        pool.map(travelCountjoinWg,options.snpfilelist)
        pool.close()
        pool.join()
        time.sleep(60)        
        filelistfile=open(options.output+'prefix.filelist','r')    
#     for SnpFile in options.snpfilelist:
#         pass
#     D_weigonSNPfilemerged=open(options.output+"merged.P123Oweigon.joinSNP","w")

#     filelistfile=open(options.output+'prefix.filelist','r')
    
    print(BBAAcountsum,ABBAcountsum,BABAcountsum)
    filelistfile.seek(0)
    filelist=filelistfile.readlines()
    for fn in filelist:
        for daftotal,n in [[delta_DerAftotal,".FreqStratifiedP1P2BBAA"],[delta_DerAftotalP1P3,'.FreqStratifiedP1P3BABA'],[delta_DerAftotalP2P3,'.FreqStratifiedP2P3ABBA']]:
            delta_DerAf=pickle.load(open(fn.strip()+n,"rb"))
            for a,b in sorted(daftotal.keys()):
                print(a,b)
                for k in sorted(daftotal[(a,b)]):
                    if k!="type1meanvalue" and k!="type2meanvalue":            
                        print(fn.strip(),len(delta_DerAf[(a,b)][k]),end="|")
                        daftotal[(a,b)][k]+=len(delta_DerAf[(a,b)][k])
                    else:
                        print(n[-4:],k,delta_DerAf[(a,b)].keys())
                        print(delta_DerAf[(a,b)][k])
                        daftotal[(a,b)][k]=numpy.mean(delta_DerAf[(a,b)][k])
    
    #BBAA
    binfile=open(options.output+"total.FreqStratifiedBBAA","w")
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(delta_DerAftotal[(a,b)]):
        print(k,end="\t",file=binfile)
    print("",file=binfile)
    for a,b in sorted(delta_DerAftotal.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAftotal[(a,b)]):
            print(delta_DerAftotal[(a,b)][k],end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()
    #BABA
    binfile=open(options.output+"total.FreqStratifiedBABA","w")
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(delta_DerAftotalP1P3[(a,b)]):
        print(k,end="\t",file=binfile)
    print("",file=binfile)
    for a,b in sorted(delta_DerAftotalP1P3.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAftotalP1P3[(a,b)]):
            print(delta_DerAftotalP1P3[(a,b)][k],end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()
    #ABBA
    binfile=open(options.output+"total.FreqStratifiedABBA","w")
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(delta_DerAftotalP2P3[(a,b)]):
        print(k,end="\t",file=binfile)
    print("",file=binfile)
    for a,b in sorted(delta_DerAftotalP2P3.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAftotalP2P3[(a,b)]):
            print(delta_DerAftotalP2P3[(a,b)][k],end="\t",file=binfile)
        print("",file=binfile)
    print(BBAAcountsum,ABBAcountsum,BABAcountsum)
#     filelistfile.seek(0)
#     fn1=filelist[0].strip()
#     rstr2=Util.random_str()
#     for fn in filelist[1:]:
#         os.system("cat "+fn1+"P123Oweigon.joinSNP "+ fn.strip()+"P123Oweigon.joinSNP > "+rstr2+"temp")
#         os.system("mv "+rstr2+"temp "+fn1+"P123Oweigon.joinSNP ")
#     print("chrNo\tpos\tANC\tDER\tP1derFreq\tP2derFreq\tP3derFreq\tP4derFreq\tBBBA\tABBA\tBABA\twigeonAF",file=D_weigonSNPfilemerged)
    binfile.close();filelistfile.close()#;D_weigonSNPfilemerged.close()
#     os.system("cat "+options.output+"merged.P123Oweigon.joinSNP "+fn1+"P123Oweigon.joinSNP ")