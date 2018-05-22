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
breaks = 5
dincrease = (maxvalue - minvalue) / breaks
delta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    delta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":[],"dxylist":[]}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
delta_DerAftotal.pop(minvalue-dincrease,minvalue);delta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":[],"dxylist":[]}

minvalue = 0
maxvalue = 1
breaks = 5
dincrease = (maxvalue - minvalue) / breaks
absdelta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    absdelta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":[],"dxylist":[]}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
absdelta_DerAftotal.pop(minvalue-dincrease,minvalue);absdelta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":[],"dxylist":[]}
ddafcaculator=Caculators.Caculate_ddaf(delta_DerAftotal,absdelta_DerAftotal)
maporder=[];infile1fstmap={}
infileref=open(options.fstwinfile,'r')
for line in infileref:
    linelist=re.split(r"\s+",line.strip())
    if  len(linelist)>5:
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
    if  len(linelist)>5:
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
            seletedregionMapByChr[regionlist[0]].append((float(regionlist[1]),float(regionlist[2])))
        else:
            seletedregionMapByChr[regionlist[0]]=[(float(regionlist[1]),float(regionlist[2]))]
    seletedtablefile.close()
    # start calculated ddaf absddaf, meanwhile, collect ddaf in high divergence and genome background
    snpfile=open(options.snpfile,'r')
    winLinAchr=[];obsexpsignalmapbychrom={};regionvalue={};backgroundvalue={}
    curchrom=re.split(r'\s+',snpfile.readline())[0]
    for line in snpfile:
        snplist=re.split(r'\s+',line.strip())
        if snplist[0]!=curchrom and curchrom!="chrNo":
            #region
            if curchrom in seletedregionMapByChr:
                ps=1;pe=1
                for s,e in seletedregionMapByChr[curchrom]:
                    win.slidWindowOverlap(winLinAchr,e,e-s,e-s,ddafcaculator,s)
                    if curchrom not in regionvalue:
                        regionvalue[curchrom]=[];backgroundvalue[curchrom]=[]
                    regionvalue[curchrom].append(copy.deepcopy(win.winValueL))
                    win.slidWindowOverlap(winLinAchr,s,s-pe,s-pe,ddafcaculator,pe+1)
                    backgroundvalue[curchrom].append(copy.deepcopy(win.winValueL))# between divergence region or before 
                    ps=s;pe=e
                else:
                    win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom]-pe,chrlenmap[curchrom]-pe,ddafcaculator,pe+1)
                    backgroundvalue[curchrom].append(copy.deepcopy(win.winValueL))# end of the last divergence to end of the chromosome. 
            else:#no divergence region
                win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom],chrlenmap[curchrom],ddafcaculator)
                backgroundvalue[curchrom]=[copy.deepcopy(win.winValueL)]
            #win
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
    for chrom in sorted(obsexpsignalmapbychrom):
        for i in range(len(obsexpsignalmapbychrom[chrom])):
            absddaf,ddaf=obsexpsignalmapbychrom[chrom][i][3]
            for a,b in sorted(absdelta_DerAftotal.keys()):#delta_DerAftotal
                if absddaf>=a and absddaf<=b:
                    fstvalue=infile1fstmap[chrom][str(i)][5]
                    dxyvalue=infile1dxymap[chrom][str(i)][5]
                    if fstvalue!="NA": absdelta_DerAftotal[(a,b)]["fstlist"].append(fstvalue)
                    if dxyvalue!="NA": absdelta_DerAftotal[(a,b)]["dxylist"].append(dxyvalue)
            for a,b in sorted(delta_DerAftotal.keys()):
                if ddaf>=a and ddaf<=b:
                    fstvalue=infile1fstmap[chrom][str(i)][5]
                    dxyvalue=infile1dxymap[chrom][str(i)][5]
                    if fstvalue!="NA": delta_DerAftotal[(a,b)]["fstlist"].append(fstvalue)
                    if dxyvalue!="NA": delta_DerAftotal[(a,b)]["dxylist"].append(dxyvalue)
    #test show result
    hdvf=open(options.output+"highdivergenceregion",'w')
    for chrom in regionvalue.keys():
        for winvaluesl in regionvalue[chrom]:
            for s,e,n,v in winvaluesl:
                print(chrom,s,e,n,v,file=hdvf)
    hdvf.close()
    bgf=open(options.output+"genomicbackgroundegion",'w')
    for chrom in backgroundvalue.keys():
        for winvaluesl in backgroundvalue[chrom]:
            for s,e,n,v in winvaluesl:
                print(chrom,s,e,n,v,file=bgf)
    bgf.close()
    binfile=open(options.output+".Stratifiedfstdxybydaf","w")
    print("bins\tbine",end="\t",file=binfile)    
    for a,b in sorted(delta_DerAftotal.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(delta_DerAftotal[(a,b)]):
            print(numpy.mean(delta_DerAftotal[(a,b)][k]),end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()

    binfile=open(options.output+".Stratifiedfstdxybyabsdaf","w")
    print("bins\tbine",end="\t",file=binfile)    
    for a,b in sorted(absdelta_DerAftotal.keys()):
        print(a,b,sep="\t",end="\t",file=binfile)
        for k in  sorted(absdelta_DerAftotal[(a,b)]):
            print(numpy.mean(absdelta_DerAftotal[(a,b)][k]),end="\t",file=binfile)
        print("",file=binfile)
    binfile.close()