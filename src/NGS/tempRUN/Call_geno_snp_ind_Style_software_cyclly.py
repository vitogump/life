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
parser.add_option("-c", "--configure", dest="configure")
parser.add_option("-t","--cmperbp",dest="cmperbp",default="3e-06")
parser.add_option("-s","--software",dest="software",help="GATK or samtools ")
parser.add_option("-1", "--ld-window-kb", dest="ldwinkb")
parser.add_option("-2", "--ld-window", dest="ldwin")
parser.add_option("-m","--sampleID_to_popmap",dest="sampleID_to_popmapfile")
parser.add_option("-d","--dilute",dest="dilute",default="1")
parser.add_option("-o","--outputpre",dest="outputpre")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

configure = open(options.configure, 'r')

sampleID_to_popmap={}
sampleID_to_popmapfile=open(options.sampleID_to_popmapfile,'r')
for line in sampleID_to_popmapfile:
    linelist=re.split(r'\s*=\s*',line.strip())
    sampleID_to_popmap[linelist[1].strip()]=linelist[0].strip()
outputprefix=options.outputpre.strip()
tempvcffile=open(outputprefix+".vcf","w")
dilute =float(options.dilute.strip())
if dilute >1 or dilute <0:
    dilute =1
chromlisttosub=configure.readlines()
print(chromlisttosub)
software=options.software.upper().strip()
cmperbp=float(options.cmperbp)

if __name__ == '__main__':
    vcfdata=VCFutil.VCF_Data(options.vcffilename.strip())
    i=0;outputfilepart=0;sumRecOfVCF=0
    for chrom in vcfdata.chromOrder:
        vcfRecOfAChrom=vcfdata.getVcfListByChrom(options.vcffilename.strip(), chrom,dilute)
        if len(vcfRecOfAChrom)<200:
            print("Call_geno_snp_ind_Style_software_cyclly","skip chrom with snps less than 100")
            continue
        else:
            sumRecOfVCF+=len(vcfRecOfAChrom)
        chrom_sub=chromlisttosub[i%len(chromlisttosub)].strip()
        if(i%len(chromlisttosub))==0 and i!=0:
            tempvcffile.close()
            VCFutil.VCF_Data.Vcf2geno_snp_ind(outputprefix+".vcf",sampleID_to_popmap,outputprefix, software,cmperbp,vcfdata.VcfIndexMap)
            tempvcffile=open(outputprefix+".vcf","w")
            break
        for pos, REF, ALT, INFO,FORMAT,samples in vcfRecOfAChrom:
            print(chrom_sub,pos,".", REF, ALT,"100",".", INFO,FORMAT,*samples,sep="\t",end="\n",file=tempvcffile)
        i+=1