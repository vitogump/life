# -*- coding: UTF-8 -*-

'''
Created on 2015-2-25

@author: liurui
'''

from optparse import OptionParser
import re, sys, os

from NGS.BasicUtil import Util, VCFutil


parser = OptionParser()
parser.add_option("-v", "--vcffile", dest="vcffilename",# action="callback",type="string",callback=useoptionvalue_previous1,
                  help="write report to FILE")
parser.add_option("-c", "--configure", dest="configure")
parser.add_option("-s","--software",dest="software",help="GATK or samtools ")
parser.add_option("-1", "--ld-window-kb", dest="ldwinkb")
parser.add_option("-2", "--ld-window", dest="ldwin")

parser.add_option("-o","--outputpre",dest="outputpre")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()

configure = open(options.configure, 'r')
cline=configure.readline()
pathtoplink = re.search(r'pathtoplink=(.*)',cline).group(1).strip()
cline=configure.readline()
temppath=re.search(r"temppath=(.*)",cline).group(1).strip()
configure.readline()
chromlisttosub=configure.readlines()
print(chromlisttosub)

software=options.software.strip()
outputprefix=options.outputpre.strip()
tempvcffile=open(outputprefix+".vcf","w");tempvcffile.close()
if __name__ == '__main__':
    vcfdata=VCFutil.VCF_Data(options.vcffilename.strip())
    if not os.path.exists(temppath):
        os.makedirs(temppath)
    os.chdir(temppath)
    i=0;outputfilepart=0
    for chrom in vcfdata.chromOrder:
        vcfRecOfAChrom=vcfdata.getVcfListByChrom(options.vcffilename.strip(), chrom)
        chrom_sub=chromlisttosub[i%len(chromlisttosub)].strip()
        if(i%len(chromlisttosub))==0:
            
            if i>0:#just for skipping the first time
                tempvcffile.close()
                VCFutil.VCF_Data.Vcf2Ped(outputprefix+".vcf",outputprefix,software,vcfdata.VcfIndexMap)
                os.system(pathtoplink+" --file "+outputprefix +" --r2 --ld-window-kb "+options.ldwinkb+" --ld-window "+options.ldwin)
                os.system("mv plink.ld plink_part"+str(outputfilepart)+".ld")
            tempvcffile=open(outputprefix+".vcf","w")
            outputfilepart+=1
        for pos, REF, ALT, INFO,FORMAT,samples in vcfRecOfAChrom:
            print(chrom_sub,pos,".", REF, ALT,"100",".", INFO,FORMAT,*samples,sep="\t",end="\n",file=tempvcffile)
        i+=1
        
               
               
    configure.close()