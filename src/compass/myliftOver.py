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

parser.add_option("-v", "--variantfilewithref", dest="variantfilewithref",action="append",nargs=4, default=[],help="vcflikefile corresponding_ref flanklen")
parser.add_option("-b", "--functionalbedlikefile", dest="functionalbedlikefile",action="append",nargs=4,default=[], help="functionalRegionfile corresponding_ref minRegionLen")
"""
functionalbedlikefile format: first line is title
chr    startpos    endpos    other info
"""
parser.add_option("-o", "--outfilename", dest="outfilename", help="it's the depth of the dir from the inputdatapath which the data file that need to be process in it,the depth of the inputdatapath is 0")
parser.add_option("-s", "--chrsignal", dest="chrsignal",default=None,help="e.g chr number follow the 'chromosome' string ")
parser.add_option("-T", "--targetREFblastdb", dest="targetREFblastdb",help="blastdb of target reference")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options,args)=parser.parse_args()
snpoutfafile=open(options.outfilename+"SNPs_flankseq.fa",'w')
bedoutfafile=open(options.outfilename+"regionsSEQ.fa",'w')
snpoutfafile.close();bedoutfafile.close()
snpoutfafile=open(options.outfilename+"SNPs_flankseq.fa",'a')
bedoutfafile=open(options.outfilename+"regionsSEQ.fa",'a')
ancestralalleletabletools=AncestralAlleletabletools(database=Util.vcfdbname, ip=Util.ip, usrname=Util.username, pw=Util.password,dbgenome=Util.genomeinfodbname)
if __name__ == '__main__':
    

    print(options.variantfilewithref)
    for chrlist,vcflikeFileName,corresponding_ref,flanklen in options.variantfilewithref:
        chromlistfile=open(chrlist,"r")
        chrmap={}
        for rec in chromlistfile:
            reclist=re.split(r'\s+',rec.strip())
            chrmap[reclist[0]]=reclist[1]        
        
        flanklen=int(flanklen)
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
            snpline=vcflinesalchr.pop(0).strip()

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
                    if chrom in chrmap:
                        chrlen=int(chrmap[chrom])
                    ancestralalleletabletools.getflankseqstooutfile(chrom, chrlen, startpostocollecteSNP, endpostocollectSNP, duckrefhandler, None, duckrefindex, flanklen, snpoutfafile, snpsOfOneChrom, None)
                    #start next chrom
                    snpsOfOneChrom=[snp]
                    chrom=snp[0]                    
                    startpostocollecteSNP=int(snp[1])
                else:#first
                    snpsOfOneChrom.append(snp);chrom=snp[0];startpostocollecteSNP=int(snp[1])
        chromlistfile.close()
    else:
        duckrefhandler.close()
        snpoutfafile.close()
            #ancestralalleletabletools.forchenli.close()
#             print()#print fa seq
    for chrlist,regionbedFName,corresponding_ref,minRegionLEN in options.functionalbedlikefile:
        chromlistfile=open(chrlist,"r")
        chrmap={}
        for rec in chromlistfile:
            reclist=re.split(r'\s+',rec.strip())
            chrmap[reclist[0]]=reclist[1]    
        minRegionLEN=int(minRegionLEN)
        regionbedf=open(regionbedFName,'r')
        duckrefhandler=open(corresponding_ref,'r')
        try:
            duckrefindex = pickle.load(open(corresponding_ref + ".myfasteridx", 'rb'))
#             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
        except IOError:
            Util.generateFasterRefIndex(corresponding_ref, corresponding_ref + ".myfasteridx",chrsignal=options.chrsignal)
            duckrefindex = pickle.load(open(corresponding_ref + ".myfasteridx", 'rb'))
        regionbedf.readline()#title
        regionsOfOneChrom=[]
        startposOfFirstREGIONs=1         
        for line in regionbedf:
            regionlist=re.split(r"\s+",line.strip())
            if chrom==regionlist[0]:
                regionsOfOneChrom.append(regionlist)
            elif regionsOfOneChrom!=[]:
                #process last
                endpostocollectREGIONs=int(regionsOfOneChrom[-1][2])
                #2，extract flank seq of variants recs
                if chrom in chrmap:
                    chrlen=int(chrmap[chrom])
                ancestralalleletabletools.getregionseqstooutfile(chrom, chrlen, startposOfFirstREGIONs, endpostocollectREGIONs, duckrefhandler, None, duckrefindex, minRegionLEN, bedoutfafile, regionsOfOneChrom, None)
                #start next chrom
                regionsOfOneChrom=[regionlist]
                chrom=regionlist[0]                    
                startposOfFirstREGIONs=int(regionlist[1])
            else:
                regionsOfOneChrom.append(regionlist);chrom=regionlist[0];startpostocollecteSNP=int(regionlist[1])
        chromlistfile.close()
    else:
        bedoutfafile.close()
        duckrefhandler.close()
    if options.variantfilewithref!=[]:
    #makeblastdb -in Setaria_italica.JGIv2.0.dna_sm.toplevel.fa -dbtype nucl -parse_seqids -out Setaria_italica.JGIv2.0.dna_sm.toplevel
        ancestralalleletabletools.callblast("blastn",options.targetREFblastdb,options.outfilename+"SNPs_flankseq.fa",options.outfilename+"SNPs_flankseq.blastout")
        ancestralalleletabletools.extarctBlastOut(options.outfilename+"SNPs_flankseq.blastout",flanklen)
    if options.functionalbedlikefile!=[] :
        ancestralalleletabletools.callblast("blastn",options.targetREFblastdb,options.outfilename+"regionsSEQ.fa",options.outfilename+"regionsSEQ.blastout") 
        ancestralalleletabletools.extarctBlastOut(options.outfilename+"regionsSEQ.blastout",flanklen)
    #ancestralalleletabletools.extarctAncestryAlleleFromBlastOut(BlastOutFile, ancestryrefFile, ancestralgenomename, ancestryrefidx, tablename, ancestralsnptable)