# -*- coding: UTF-8 -*-
import re,sys
from optparse import OptionParser
from NGS.BasicUtil import Util, VCFutil
'''
Created on 2015-4-22

@author: liurui
'''
parser = OptionParser()
parser.add_option("-v", "--vcffile", dest="vcffilename",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--configure", dest="configure",default=None,help="not neccery")
parser.add_option("-l", "--chromlistfilename", dest="chromlistfilename", help="i")

parser.add_option("-t","--Morganperbp",dest="Morganperbp",default="5.4696217209617786e-8")
parser.add_option("-s","--software",dest="software",help="GATK or samtools ")
parser.add_option("-1", "--ld-window-kb", dest="ldwinkb")
parser.add_option("-2", "--ld-window", dest="ldwin")
parser.add_option("-m","--sampleID_to_popmap",dest="sampleID_to_popmapfile")
parser.add_option("-d","--dilute",dest="dilute",default="1")
parser.add_option("-D","--dilutetodensity",dest="dilutetodensity",default=None,help="snp per kb")
parser.add_option("-o","--outputpre",dest="outputpre")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()
if options.configure!=None:
    configure = open(options.configure, 'r')
    chromlisttosub=configure.readlines()
else:
    chromlisttosub=None

sampleID_to_popmap={}

sampleID_to_popmapfile=open(options.sampleID_to_popmapfile,'r')
for line in sampleID_to_popmapfile:
    linelist=re.split(r'\s*=\s*',line.strip())
    sampleID_to_popmap[linelist[1].strip()]=linelist[0].strip()
outputprefix=options.outputpre.strip()

if options.dilutetodensity!=None and options.dilute=="1":
    dilutetodensity=float(options.dilutetodensity.strip())
    dilute =1
elif options.dilute!="1" and options.dilutetodensity==None:
    dilute =float(options.dilute.strip())
    if dilute >1 or dilute <0:
        dilute =1
else:
    print("error")
    exit(-1)

print(chromlisttosub)
software=options.software.upper().strip()
Morganperbp=float(options.Morganperbp)
chromlistfile=open(options.chromlistfilename,"r")
chromlist=[]
for chrrow in chromlistfile:
    chrrowlist=re.split(r'\s+',chrrow.strip())
    chromlist.append(chrrowlist[0].strip())
tempvcffile=open(outputprefix+".vcf","w")
if __name__ == '__main__':
    vcfdata=VCFutil.VCF_Data(options.vcffilename.strip())
    i=0;outputfilepart=0;sumRecOfVCF=0
    if chromlisttosub==None:
        lastpos=0
    for chrom in chromlist:
        if chrom not in vcfdata.chromOrder:
            continue
        vcfRecOfAChrom=vcfdata.getVcfListByChrom(chrom,dilute,dilutetodensity=dilutetodensity)
        if len(vcfRecOfAChrom)<30:
            print("Call_geno_snp_ind_Style_software_cyclly","skip chrom with snps less than 100")
            continue
        else:
            sumRecOfVCF+=len(vcfRecOfAChrom)
        if chromlisttosub!=None:
            chrom_sub=chromlisttosub[i%len(chromlisttosub)].strip()
        else:
            chrom_sub="1"#chrom
        if chromlisttosub!=None and (i%len(chromlisttosub))==0 and i!=0:
            tempvcffile.close()
            VCFutil.VCF_Data.Vcf2geno_snp_ind(outputprefix+".vcf",sampleID_to_popmap,outputprefix, software,Morganperbp,"all",vcfdata.VcfIndexMap)
            tempvcffile=open(outputprefix+".vcf","w")
            outputfilepart+=1
        for pos, REF, ALT, INFO,FORMAT,samples in vcfRecOfAChrom:
            if chromlisttosub==None:
                pos=lastpos+pos
            print(chrom_sub,pos,".", REF, ALT,"100",".", INFO,FORMAT,*samples,sep="\t",end="\n",file=tempvcffile)
        else:
            lastpos=pos
        i+=1
    else:
        if chromlisttosub==None:
            tempvcffile.close()
            VCFutil.VCF_Data.Vcf2geno_snp_ind(outputprefix+".vcf",sampleID_to_popmap,outputprefix, software,Morganperbp,"all",vcfdata.VcfIndexMap)
    print("finish")