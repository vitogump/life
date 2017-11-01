# -*- coding: UTF-8 -*-
'''
Created on 2017年10月24日

@author: liurui
'''
from optparse import OptionParser
import pickle, re, sys

from NGS.BasicUtil import Util
from src.NGS.Service.Ancestralallele import AncestralAlleletabletools


parser = OptionParser()

#"output data name is defined as 'inputdatapath folder name'+'is subfolder name'+'is subfolder name'+..."

parser.add_option("-v", "--variantfilewithref", dest="variantfilewithref",action="append",nargs=2, help="vcflikefile corresponding_ref")
parser.add_option("-r", "--objectref", dest="objectref", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-s", "--chrsignal", dest="chrsignal",help="chromosome")
parser.add_option("-f", "--flanklen", dest="flanklen",default='70',help="ref fa file mode2")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
outfile=open('test.txt','w')
flanklen=int(options.flanklen)
ancestralalleletabletools=AncestralAlleletabletools(database=Util.vcfdbname, ip=Util.ip, usrname=Util.username, pw=Util.password,dbgenome=Util.genomeinfodbname)
if __name__ == '__main__':
    
    
    for vcflikeFileName,corresponding_ref in options.variantfilewithref:
        duckrefhandler=open(corresponding_ref,'r')
        try:
            duckrefindex = pickle.load(open(corresponding_ref + ".myfasteridx", 'rb'))
#             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        except IOError:
            Util.generateFasterRefIndex(corresponding_ref, corresponding_ref + ".myfasteridx",chrsignal=options.chrsignal)
            duckrefindex = pickle.load(open(corresponding_ref + ".myfasteridx", 'rb'))
        vcflikefile=open(vcflikeFileName,'r')
        vcflinesalchr=vcflikefile.readlines()
        #1，read variations
        chrom=None
        snpsOfOneChrom=[]
        startpostocollecteSNP=1        
        while vcflinesalchr:
            snpline=vcflinesalchr.pop(0)
            print(snpline)
            if snpline[0]=="#" or snpline.lower().find("chrom")==0:#title
                continue
            else:
                snp=re.split(r"\s+",snpline.strip())
                
                if chrom==snp[0]:
                    snpsOfOneChrom.append(snp)
                    
                elif snpsOfOneChrom!=[]:
                    #process last

                    endpostocollectSNP=int(snpsOfOneChrom[-1][1])
                    #2，extract flank seq of variants recs
                    ancestralalleletabletools.getflankseqs(chrom, None, startpostocollecteSNP, endpostocollectSNP, duckrefhandler, None, duckrefindex, flanklen, outfile, snpsOfOneChrom, None)
                    #start next chrom
                    snpsOfOneChrom=[snp]
                    chrom=snp[0]                    
                    startpostocollecteSNP=int(snp[1])
                else:#first
                    snpsOfOneChrom.append(snp);chrom=snp[0];startpostocollecteSNP=int(snp[1])
            
        else:
            duckrefhandler.close()
            outfile.close()
                
            print()#print fa seq
#     ancestralalleletabletools.callblat()
#     ancestralalleletabletools.extarctAncestryAlleleFromBlastOut(BlastOutFile, ancestryrefFile, ancestralgenomename, ancestryrefidx, tablename, ancestralsnptable)