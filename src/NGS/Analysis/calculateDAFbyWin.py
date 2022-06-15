'''
Created on 2018-5-20

@author: Dr.liu
'''
from src.NGS.BasicUtil import Caculators
from src.NGS.BasicUtil.Util import Window
from optparse import OptionParser
import re,copy,numpy


parser = OptionParser()
parser.add_option("-s", "--snpfile", dest="snpfile", help="early,pairfst,pbs,lsbl,is")
parser.add_option("-b", "--seletedtable", dest="seletedtable")
parser.add_option("-c", "--chrfile", dest="chrfile")
parser.add_option("-f", "--fstwinfile", dest="fstwinfile")
parser.add_option("-d", "--dxywinfile", dest="dxywinfile")
parser.add_option("-o", "--output", dest="output")

(options, args) = parser.parse_args()

print(options.output)
chrlenmap={}
f=open(options.chrfile,'r')
for line in f:
    linelist=re.split(r'\s+',line.strip())
    chrlenmap[linelist[0]]=int(linelist[1])
f.close()

win=Window()

minvalue = -1
maxvalue = 1
breaks = 10
dincrease = (maxvalue - minvalue) / breaks
delta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    delta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":[],"dxylist":[],"snpcount":[]}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
delta_DerAftotal.pop(minvalue-dincrease,minvalue);delta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":[],"dxylist":[],"snpcount":[]}

minvalue = 0
maxvalue = 1
breaks = 5
dincrease = (maxvalue - minvalue) / breaks
absdelta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    absdelta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":[],"dxylist":[],"snpcount":[]}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
absdelta_DerAftotal.pop(minvalue-dincrease,minvalue);absdelta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":[],"dxylist":[],"snpcount":[]}
ddafcaculator=Caculators.Caculate_ddaf(delta_DerAftotal,absdelta_DerAftotal)
maporder=[];infile1fstmap={}
infileref=open(options.fstwinfile,'r')
for line in infileref:
    linelist=re.split(r"\s+",line.strip())
    if  len(linelist)>5 and False:
        if linelist[5]=="na" or linelist[5]=="NA" or re.search(r"inf", linelist[5])!=None or linelist[6]=="na" or linelist[6]=="NA" or re.search(r"inf", linelist[5])!=None:
            print("skip",linelist)
            continue
    maporder.append((linelist[0].strip(),linelist[1].strip()))
    try:
        infile1fstmap[linelist[0].strip()][linelist[1].strip()]=linelist
    except KeyError:
        infile1fstmap[linelist[0].strip()]={linelist[1].strip():linelist}
infileref.close()  
infileref=open(options.dxywinfile,'r');infile1dxymap={}
for line in infileref:
    linelist=re.split(r"\s+",line.strip())
    if  len(linelist)>5 and False:
        if linelist[5]=="na" or linelist[5]=="NA" or re.search(r"inf", linelist[5])!=None or linelist[6]=="na" or linelist[6]=="NA" or re.search(r"inf", linelist[5])!=None:
            print("skip",linelist)
            continue
#     maporder.append((linelist[0].strip(),linelist[1].strip()))
    try:
        infile1dxymap[linelist[0].strip()][linelist[1].strip()]=linelist
    except KeyError:
        infile1dxymap[linelist[0].strip()]={linelist[1].strip():linelist}
infileref.close() 
if __name__ == '__main__':
    #read high heterogeneous region
    seletedtablefile=open(options.seletedtable,'r')
    print(seletedtablefile.readline())
    seletedregionMapByChr={}
    for line in seletedtablefile:
        regionlist=re.split(r'\s+',line.strip())
        if regionlist[0] in seletedregionMapByChr:
            seletedregionMapByChr[regionlist[0]].append((int(regionlist[1]),int(regionlist[2])))
        else:
            seletedregionMapByChr[regionlist[0]]=[(int(regionlist[1]),int(regionlist[2]))]
    seletedtablefile.close()
    print(seletedregionMapByChr)
    # start calculated ddaf absddaf, meanwhile, collect ddaf in high divergence and genome background
    snpfile=open(options.snpfile,'r')
    winLinAchr=[];obsexpsignalmapbychrom={};regionvalue={};backgroundvalue={}
    curchrom=re.split(r'\s+',snpfile.readline())[0]
    for line in snpfile:
        snplist=re.split(r'\s+',line.strip())
        if snplist[0] not in chrlenmap:continue
        if snplist[0]!=curchrom and curchrom!="chrNo":
            #region
            backgroundvalue[curchrom]=[]
            if curchrom in seletedregionMapByChr:
                ps=1;pe=1;regionvalue[curchrom]=[]
                for s,e in seletedregionMapByChr[curchrom]:
                    if (winLinAchr[0][0]>s and winLinAchr[0][0]<e) or (winLinAchr[-1][0]>s and winLinAchr[-1][0]<e) or (winLinAchr[0][0]<=s and winLinAchr[-1][0]>=e):
                        print(curchrom,s,e,file=open("adsfasdf",'a'))
                        print(curchrom)
                        win.slidWindowOverlap(winLinAchr,e,e-s+1,e-s+1,ddafcaculator,s)
                        regionvalue[curchrom].append(copy.deepcopy(win.winValueL))
                    if s>pe+1 and winLinAchr[0][0]>pe and winLinAchr[-1][0]<s:
                        win.slidWindowOverlap(winLinAchr,s,s-pe+1,s-pe+1,ddafcaculator,pe)
                        backgroundvalue[curchrom].append(copy.deepcopy(win.winValueL))# between divergence region or before 
                    ps=s;pe=e
                else:
                    if chrlenmap[curchrom]<=pe+1 or winLinAchr[-1][0]<=pe:
                        print(curchrom,chrlenmap[curchrom],"or",winLinAchr[-1][0],"<",pe+1)
                    else:
                        win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom]-pe+1,chrlenmap[curchrom]-pe+1,ddafcaculator,pe)
                        backgroundvalue[curchrom].append(copy.deepcopy(win.winValueL))# end of the last divergence to end of the chromosome. 
            else:#no divergence region
                win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom]+1,chrlenmap[curchrom]+1,ddafcaculator)
                backgroundvalue[curchrom]=[copy.deepcopy(win.winValueL)]
            #win
#             print(len(winLinAchr),winLinAchr[0],curchrom)
            win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],40000,20000,ddafcaculator)
            obsexpsignalmapbychrom[curchrom]=copy.deepcopy(win.winValueL)
            curchrom=snplist[0];snplist[1]=int(snplist[1])
            winLinAchr=[snplist[1:]]
        elif curchrom!="chrNo":
            snplist[1]=int(snplist[1]);winLinAchr.append(snplist[1:])
        else:#first time
            curchrom=snplist[0];snplist[1]=int(snplist[1])
            winLinAchr=[snplist[1:]]
    # stratified fst,dxy by ddaf bin
    f=open("correlation.fstdxy_daf",'w')
    for chrom in sorted(obsexpsignalmapbychrom):
        for i in range(len(obsexpsignalmapbychrom[chrom])):
            absddaf,ddaf=obsexpsignalmapbychrom[chrom][i][3];sitscountforddaf=obsexpsignalmapbychrom[chrom][i][2]
            print(infile1dxymap[chrom][str(i)][5])
            if sitscountforddaf!=0: print(chrom,i,infile1fstmap[chrom][str(i)][5],infile1dxymap[chrom][str(i)][5],absddaf/sitscountforddaf,ddaf/sitscountforddaf,file=f);meanabsddaf=absddaf/sitscountforddaf;meanddaf=ddaf/sitscountforddaf
            for a,b in sorted(absdelta_DerAftotal.keys()):#delta_DerAftotal
                if meanabsddaf>=a and meanabsddaf<=b:
                    fstvalue=infile1fstmap[chrom][str(i)][5]
                    dxyvalue=infile1dxymap[chrom][str(i)][5]
                    absdelta_DerAftotal[(a,b)]["snpcount"].append(sitscountforddaf)
                    if fstvalue!="NA": absdelta_DerAftotal[(a,b)]["fstlist"].append(float(fstvalue))
                    if dxyvalue!="NA": absdelta_DerAftotal[(a,b)]["dxylist"].append(float(dxyvalue))
            for a,b in sorted(delta_DerAftotal.keys()):
                if meanddaf>=a and meanddaf<=b:
                    fstvalue=infile1fstmap[chrom][str(i)][5]
                    dxyvalue=infile1dxymap[chrom][str(i)][5]
                    delta_DerAftotal[(a,b)]["snpcount"].append(sitscountforddaf)
                    if fstvalue!="NA": delta_DerAftotal[(a,b)]["fstlist"].append(float(fstvalue))
                    if dxyvalue!="NA": delta_DerAftotal[(a,b)]["dxylist"].append(float(dxyvalue))
    f.close()
    #test show result
    hdvf=open(options.output+"highdivergenceregion",'w')
    print("chrNo\tstartpos\tendpos\tabsddaf\tddaf",file=hdvf)
    for chrom in regionvalue.keys():
        for winvaluesl in regionvalue[chrom]:
            for s,e,n,v in winvaluesl:
                print(chrom,s,e,n,*v,sep="\t",file=hdvf)
    hdvf.close()
    
    bgf=open(options.output+"genomicbackgroundegion",'w')
    print("chrNo\tstartpos\tendpos\tabsddaf\tddaf",file=bgf)
    for chrom in backgroundvalue.keys():
        for winvaluesl in backgroundvalue[chrom]:
            for s,e,n,v in winvaluesl:
                print(chrom,s,e,n,*v,sep="\t",file=bgf)
    bgf.close()
    binfile=open(options.output+".Stratifiedfstdxybydaf","w")
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(delta_DerAftotal[(a,b)]):
        
        print(k,end="\t",file=binfile)
    print("",file=binfile) 
    for a,b in sorted(delta_DerAftotal.keys()):
        print(*delta_DerAftotal[(a,b)]["dxylist"],sep="\n",file=open("dxyddaf.winvalue",'a'))
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAftotal[(a,b)]):
            print(numpy.mean(delta_DerAftotal[(a,b)][k]),end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()

    binfile=open(options.output+".Stratifiedfstdxybyabsdaf","w")
    print("bins\tbine",end="\t",file=binfile)
    for k in sorted(absdelta_DerAftotal[list(absdelta_DerAftotal.keys())[0]]):
        
        print(k,end="\t",file=binfile)
    print("",file=binfile)
    for a,b in sorted(absdelta_DerAftotal.keys()):
        print(*absdelta_DerAftotal[(a,b)]["dxylist"],sep="\n",file=open("dxyabsddaf.winvalue",'a'))
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(absdelta_DerAftotal[(a,b)]):
            print(numpy.mean(absdelta_DerAftotal[(a,b)][k]),end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()