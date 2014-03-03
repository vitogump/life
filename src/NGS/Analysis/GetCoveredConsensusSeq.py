# -*- coding: UTF-8 -*-
from NGS.BasicUtil import *
from optparse import OptionParser
import pickle
import re
import sys

'''
Created on 2014-2-12

@author: liurui
'''
mindepth=2
#if len(sys.argv) != 6:
#    print(len(sys.argv))
#    print("python GetConveredConsensusSeq.py [ref.fa] [gtffile] [vcffile(withheader)] [genomeCoveragefile] [outfileprename] ")
#    exit(-1)
parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa",
                  help="reference.fa", metavar="FILE")
parser.add_option("-g", "--gtffile", dest="gtffile", help="gtffile")
parser.add_option("-v", "--vcffile", dest="variants", help="variants")
parser.add_option("-d", "--genomeCoveragefile", dest="genomedepth", help="genomeCoveragefile")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-C", "--whetherwriteconsensus", dest="cnsornot", action='store_false', help="use this option means you will not print the cns.fa file")
parser.add_option("-n", "--speciesname", dest="species", help="species name")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
reffa = open(options.reffa, 'r')
#gtffile = open(options.gtffile, 'r')
vcffile = open(options.variants, 'r')
covfile = open(options.genomedepth, 'r')


cns_string=""
aa_string=""
cdscns_string=""
outcns = open(options.outfileprename + "_cns.fa", 'w')
outaa = open(options.outfileprename + "_aa.fa", 'w')
outcdscns = open(options.outfileprename + "_cdscns.fa", 'w')
cdsmap={}
if __name__ == '__main__':
    depthfile = Util.GATK_depthfile(options.genomedepth, options.genomedepth + ".index")
    species_idx = depthfile.title.index("Depth_for_" + options.species)
    vcfpop = VCFutil.VCF_Data(options.variants)  # new a class
    RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa)
    gtfMap=Util.getGtfMap(options.gtffile)
    
    lastposofdepthfp = 0#because this time RefSeqMap[0] is 0
    while currentChromNO != "end of the reffile":
        currentBaselocinGenome = RefSeqMap[currentChromNO][0]+1
        statue = depthfile.set_depthfilefp(currentChromNO, currentBaselocinGenome, lastposofdepthfp)
        nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
        frontmostpos=nearestGenes.geneOverlapList[0][2];Rearmostpos=nearestGenes.geneOverlapList[-1][3]
        if len(RefSeqMap[currentChromNO])<=(Rearmostpos-frontmostpos+1):
            curposoffilehandler=reffa.tell()
            reffa_suplemtry= open(options.reffa, 'r')
            reffa_suplemtry.seek(curposoffilehandler)
            RefSeqMap_suplemtry, lastchromNo_suplemtry = Util.getRefSeqMap(reffa_suplemtry, currentChromNO=currentChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
            RefSeqMap[currentChromNO]+=RefSeqMap_suplemtry[1:]
            reffa_suplemtry.close()
            nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
        #the the use of if block upside is that make sure RefSeqMap[currentChromNO] has enough seq contain the geneOverlapList scope
        vcflist_A_chrom = vcfpop.getVcfListByChrom(options.variants, currentChromNO)
        idx_vcf=0
        idx_RefSeq=1
        while idx_RefSeq!=len(RefSeqMap[currentChromNO]):
            if currentBaselocinGenome<frontmostpos:
                depth_chrom, depth_pos, depth_line, depth_linelist = depthfile.getnextposline()    
                if depth_chrom != currentChromNO or depth_pos!=currentBaselocinGenome:
                    statue = depthfile.set_depthfilefp(currentChromNO, currentBaselocinGenome, lastposofdepthfp)
                    if statue =="didn't find":
                        print("warning:"+options.genomedepth+"didn't have this genome postion:"+currentChromNO+" , "+currentBaselocinGenome)
                        cns_string+="n"
                        idx_RefSeq+=1
                        continue
#                else:# normal situation. usually
                if int(depth_linelist[species_idx]) >= mindepth:
                    if vcflist_A_chrom[idx_vcf][0]==currentBaselocinGenome:
                        cns_string+=vcflist_A_chrom[idx_vcf][2]
                        idx_RefSeq+=(len(vcflist_A_chrom[idx_vcf][1])-1)
                        idx_vcf+=1
                    else:
                        cns_string+=RefSeqMap[currentChromNO][idx_RefSeq]
                        idx_RefSeq+=1
                else:
                    cns_string+="n"
                    idx_RefSeq+=1
                        
            else:
                cds_map,aa_map,cns_append,idx_vcf=nearestGenes.getgeneConsensus(RefSeqMap[currentChromNO], idx_RefSeq, vcflist_A_chrom, idx_vcf, depthfile)
                cns_string+=cns_append
                for geneName in cds_map.keys():#write to file
                    print(">"+geneName+"\n",file=outcdscns)
                    print(">"+geneName+"\n",file=outaa)
                    i =0#write cds seq to file
                    cdsstrline=cds_map[geneName][i:i+60]
                    while len(cdsstrline)==60:
                        print(cdsstrline+"\n",file=outcdscns)
                        i+=60
                        cdsstrline=cds_map[geneName][i:i+60]
                    else:
                        print(cdsstrline+"\n",file=outcdscns)
                    i =0#write protein seq to file
                    aastrline=aa_map[geneName][i:i+60]
                    while len(aastrline)==60:
                        print(aastrline+"\n",file=outaa)
                        i+=60
                        aastrline=aa_map[geneName][i:i+60]
                    else:
                        print(aastrline+"\n",file=outaa)    
                idx_RefSeq+=(Rearmostpos-frontmostpos+1);currentBaselocinGenome+=(Rearmostpos-frontmostpos+1)
                nearestGenes=Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
                frontmostpos=nearestGenes.geneOverlapList[0][2];Rearmostpos=nearestGenes.geneOverlapList[-1][3]
                if len(RefSeqMap[currentChromNO])<=(Rearmostpos-frontmostpos+1):
                    curposoffilehandler=reffa.tell()
                    reffa_suplemtry= open(options.reffa, 'r')
                    reffa_suplemtry.seek(curposoffilehandler)
                    RefSeqMap_suplemtry, lastchromNo_suplemtry = Util.getRefSeqMap(reffa_suplemtry, currentChromNO=currentChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
                    RefSeqMap[currentChromNO]+=RefSeqMap_suplemtry[1:]
                    reffa_suplemtry.close()
                    nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
        else:
            lastposofdepthfp = depthfile.depthfilefp.tell()
            print(cns_string,end="",file=outcns)
            if nextChromNO == currentChromNO:
                cns_string=""
            else:
                cns_string="\n>"+nextChromNO+"\n"
            RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=nextChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
    else:
        print("finish")        
        
        
#                    if not nearestGenes.geneOverlapList:#empty
#                        nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
#                        frontmostpos=nearestGenes.geneOverlapList[0][2];Rearmostpos=nearestGenes.geneOverlapList[-1][3]
#                        if len(RefSeqMap[currentChromNO])<=(Rearmostpos-frontmostpos+1):
#                            curposoffilehandler=reffa.tell()
#                            reffa_suplemtry= open(options.reffa, 'r')
#                            reffa_suplemtry.seek(curposoffilehandler)
#                            RefSeqMap_suplemtry, lastchromNo_suplemtry = Util.getRefSeqMap(reffa_suplemtry, currentChromNO=currentChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
#                            RefSeqMap[currentChromNO]+=RefSeqMap_suplemtry[1:]
#                            reffa_suplemtry.close()
#                            nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome,RefSeqMap[currentChromNO])
#                        #the the use of if block upside is that make sure RefSeqMap[currentChromNO] has enough seq contain the geneOverlapList scope
#                    for gene in nearestGenes.geneOverlapList:
#                        if currentBaselocinGenome>=gene[1] and currentBaselocinGenome<=gene[2]:
#                        Util.getCDS()
#                else:
#                    cns_string+="n"
#            idx_RefSeq+=1
#
#        lastposofdepthfp = depthfile.depthfilefp.tell()
#        print(cns_string,end="",file=outcns)
#        if nextChromNO == currentChromNO:
#            cns_string=""
#        else:
#            cns_string="/n>"+nextChromNO+"/n"
#        RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=nextChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
#    else:
#        pass
            
    
#gtffile.close()    
reffa.close()
vcffile.close()
covfile.close()    
outcns.close()
outaa.close()
outcdscns.close()

