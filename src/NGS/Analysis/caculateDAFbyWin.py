'''
Created on 2018-5-20

@author: Dr.liu
'''
from src.NGS.BasicUtil import Caculators
from src.NGS.BasicUtil.Util import Window
from optparse import OptionParser
import re,copy


parser = OptionParser()
parser.add_option("-s", "--snpfile", dest="snpfile", help="early,pairfst,pbs,lsbl,is")
parser.add_option("-b", "--seletedtable", dest="seletedtable")
parser.add_option("-c", "--chrfile", dest="chrfile")
parser.add_option("-f", "--fstwinfile", dest="fstwinfile")
parser.add_option("-d", "--dxywinfile", dest="dxywinfile")
(options, args) = parser.parse_args()
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
    delta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":0,"dxylist":0}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
delta_DerAftotal.pop(minvalue-dincrease,minvalue);delta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":0,"dxylist":0}

minvalue = 0
maxvalue = 1
breaks = 5
dincrease = (maxvalue - minvalue) / breaks
absdelta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    absdelta_DerAftotal[(minvalue,minvalue + dincrease)]={"fstlist":0,"dxylist":0}
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
absdelta_DerAftotal.pop(minvalue-dincrease,minvalue);absdelta_DerAftotal[(minvalue-dincrease,maxvalue)]={"fstlist":0,"dxylist":0}
ddafcaculator=Caculators.Caculate_ddaf(delta_DerAftotal,absdelta_DerAftotal)
maporder=[];infile1fstmap={}
infileref=open(options.fstwinfile,'r')
for line in infileref:
    linelist=re.split(r"\s+",line.strip())
    if options.rmna and len(linelist)>5:
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
    if options.rmna and len(linelist)>5:
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
    # start calculated ddaf 
    snpfile=open(options.snpfile,'r')
    winLinAchr=[];obsexpsignalmapbychrom={};regionvalue={};backgroundvalue={}
    curchrom=re.split(r'\s+',snpfile.readline())[0]
    for line in snpfile:
        snplist=re.split(r'\s+',line.strip())
        if snplist[0]!=curchrom and curchrom!="chrNo":
            #region
            if curchrom in seletedregionMapByChr:
                ps=0;pe=0
                for s,e in seletedregionMapByChr[curchrom]:
                    win.slidWindowOverlap(winLinAchr,e,e-s,e-s,ddafcaculator,s)
                    regionvalue[curchrom]=copy.deepcopy(win.winValueL)
                    win.slidWindowOverlap(winLinAchr,s,s-pe,s-pe,ddafcaculator,pe)
                    backgroundvalue[curchrom]=copy.deepcopy(win.winValueL)# between divergence region or before 
                    ps=s;pe=e
                else:
                    win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom]-pe,chrlenmap[curchrom]-pe,ddafcaculator,pe)
                    backgroundvalue[curchrom]=copy.deepcopy(win.winValueL)# end of the last divergence to end of the chromosome. 
            else:#no divergence region
                win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],chrlenmap[curchrom],chrlenmap[curchrom],ddafcaculator)
                backgroundvalue[curchrom]=copy.deepcopy(win.winValueL)
            #win
            win.slidWindowOverlap(winLinAchr,chrlenmap[curchrom],40000,20000,ddafcaculator)
            obsexpsignalmapbychrom[curchrom]=copy.deepcopy(win.winValueL)
            curchrom=snplist[0]
            winLinAchr=[snplist[1:]]
        elif curchrom!="chrNo":
            winLinAchr.append(snplist[1:])
        else:#first time
            curchrom=snplist[0]
            winLinAchr=[snplist[1:]]
    # stratified by ddaf bin
    for chrom in sorted(obsexpsignalmapbychrom):
        for i in range(len(obsexpsignalmapbychrom[chrom])):
            absddaf,ddaf=obsexpsignalmapbychrom[chrom][i][3]
            for a,b in sorted(absdelta_DerAftotal.keys()):#delta_DerAftotal
                if absddaf>a and absddaf<=b:
    for chrom in obsexpsignalmapbychrom.keys():
        for s,e,n,v in backgroundvalue
        print(chrom+"\t"+)
        