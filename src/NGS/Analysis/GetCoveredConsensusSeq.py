# -*- coding: UTF-8 -*-

from NGS.BasicUtil import Util, VCFutil
from optparse import OptionParser
import re
import sys,pickle

'''
Created on 2014-2-12

@author: liurui
'''
mindepth = 2
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

with open(options.reffa,'r') as testminintervalbetweengenes_basesperfaline:
    testminintervalbetweengenes_basesperfaline.readline()
    minintervalbetweengenes_basesperfaline=len(testminintervalbetweengenes_basesperfaline.readline().strip())
    print(minintervalbetweengenes_basesperfaline)
#gtffile = open(options.gtffile, 'r')
vcffile = open(options.variants, 'r')
#covfile = open(options.genomedepth, 'r')


cns_string = ">"
aa_string = ""
cdscns_string = ""
outcns = open(options.outfileprename + "_cns.fa", 'w')
outaa = open(options.outfileprename + "_aa.fa", 'w')
outcdscns = open(options.outfileprename + "_cdscns.fa", 'w')
cdsmap = {}
if __name__ == '__main__':
    if options.genomedepth !=None:
        depthfile = Util.GATK_depthfile(options.genomedepth, options.genomedepth + ".index")
        species_idx = depthfile.title.index("Depth_for_" + options.species)
        Considerdepth=True
    else:
        Considerdepth=False
        depthfile=None
        species_idx=-1
    vcfpop = VCFutil.VCF_Data(options.variants)  # new a class
    RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa)
    print(currentChromNO,nextChromNO)
    cns_string+=currentChromNO+"\n"
    gtfMap = Util.getGtfMap(options.gtffile)
    
    lastposofdepthfp = 0#because this time RefSeqMap[0] is 0
    vcfchrom="begin"
    while currentChromNO != "end of the reffile":
        print("\t\twhile loop:", currentChromNO)
        currentBaselocinGenome = RefSeqMap[currentChromNO][0] + 1
#        statue = depthfile.set_depthfilefp(currentChromNO, currentBaselocinGenome, lastposofdepthfp)
#        depth_chrom, depth_pos, depth_linelist,lastposofdepthfp = depthfile.getnextposline()
        if currentChromNO in gtfMap:
            gtfListOfCurrentChrom = gtfMap[currentChromNO]
        else:
            gtfListOfCurrentChrom = None
        nearestGenes = Util.genes(gtfListOfCurrentChrom, currentBaselocinGenome, RefSeqMap[currentChromNO])
        if nearestGenes.geneOverlapList:
            frontmostpos = nearestGenes.geneOverlapList[0][2];Rearmostpos = nearestGenes.geneOverlapList[-1][3]
        if nearestGenes.geneOverlapList and len(RefSeqMap[currentChromNO])-1+RefSeqMap[currentChromNO][0]+1 < Rearmostpos:
            bases_to_get=Rearmostpos -len(RefSeqMap[currentChromNO])-RefSeqMap[currentChromNO][0]
            lins_to_get=int(bases_to_get/minintervalbetweengenes_basesperfaline)+1
            curposoffilehandler = reffa.tell()
            reffa_suplemtry = open(options.reffa, 'r')
            reffa_suplemtry.seek(curposoffilehandler)
            RefSeqMap_suplemtry, currentchromNo_suplemtry,nextChromNO = Util.getRefSeqMap(reffa_suplemtry, currentChromNO=currentChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1,linesOnce=lins_to_get)
            if currentchromNo_suplemtry!=currentChromNO:
                print("this if block just for test")
                exit(-1)
            RefSeqMap[currentChromNO] += RefSeqMap_suplemtry[currentchromNo_suplemtry][1:]
            reffa_suplemtry.close()
            nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome, RefSeqMap[currentChromNO])
        #the the use of if block upside is that make sure RefSeqMap[currentChromNO] has enough seq contain the geneOverlapList scope
        if vcfchrom!=currentChromNO:
            vcflist_A_chrom = vcfpop.getVcfListByChrom(currentChromNO)
            vcfchrom=currentChromNO
            if vcflist_A_chrom:
                idx_vcf=0
            else:
                idx_vcf=-1

        idx_RefSeq = 1
        while idx_RefSeq < len(RefSeqMap[currentChromNO]):
            if (not nearestGenes.geneOverlapList) or currentBaselocinGenome < frontmostpos:
                if Considerdepth:
                    depth_linelist = depthfile.getdepthByPos(currentChromNO, currentBaselocinGenome)

                if (not Considerdepth) or( depth_linelist and int(depth_linelist[species_idx]) >= mindepth):
                    if vcflist_A_chrom and idx_vcf < len(vcflist_A_chrom) and vcflist_A_chrom[idx_vcf][0] == currentBaselocinGenome:
                        if re.search(r'[^a-zA-Z]', vcflist_A_chrom[idx_vcf][2]) != None:#contain ',' ie. multiple alle
                            cns_string += vcflist_A_chrom[idx_vcf][1]
                            idx_RefSeq += len(vcflist_A_chrom[idx_vcf][1]);currentBaselocinGenome += len(vcflist_A_chrom[idx_vcf][1])
                            idx_vcf += 1 
                        else:
                            cns_string += vcflist_A_chrom[idx_vcf][2]
                            idx_RefSeq += len(vcflist_A_chrom[idx_vcf][1]);currentBaselocinGenome += len(vcflist_A_chrom[idx_vcf][1]) 
                            idx_vcf += 1
                    else:
                        cns_string += RefSeqMap[currentChromNO][idx_RefSeq]
                        idx_RefSeq += 1;currentBaselocinGenome += 1
                else:
                    cns_string += "N"
                    idx_RefSeq += 1;currentBaselocinGenome += 1
                        
            else:
                print(idx_RefSeq,currentBaselocinGenome)
                cds_map, aa_map, cns_append, idx_vcf = nearestGenes.getgeneConsensus(RefSeqMap[currentChromNO], idx_RefSeq, vcflist_A_chrom, idx_vcf, depthfile)
                cns_string += cns_append
                for geneName in cds_map.keys():#write to file
                    print(">" + geneName ,end="\n", file=outcdscns)
                    print(">transcript:" + geneName ,end="\n", file=outaa)
                    i = 0#write cds seq to file
                    cdsstrline = "".join(cds_map[geneName][i:i + 60])
                    while len(cdsstrline) == 60:
                        print(cdsstrline ,end="\n", file=outcdscns)
                        i += 60
                        cdsstrline = "".join(cds_map[geneName][i:i + 60])
                    else:
                        print(cdsstrline ,end="\n", file=outcdscns)
                    i = 0#write protein seq to file
                    aastrline = "".join(aa_map[geneName][i:i + 60])
                    while len(aastrline) == 60:
                        print(aastrline ,end="\n", file=outaa)
                        i += 60
                        aastrline = "".join(aa_map[geneName][i:i + 60])
                    else:
                        print(aastrline ,end="\n", file=outaa)    
                idx_RefSeq += (Rearmostpos - frontmostpos + 1);currentBaselocinGenome += (Rearmostpos - frontmostpos + 1)
                nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome, RefSeqMap[currentChromNO])
                if nearestGenes.geneOverlapList:
                    frontmostpos = nearestGenes.geneOverlapList[0][2];Rearmostpos = nearestGenes.geneOverlapList[-1][3]
                if nearestGenes.geneOverlapList and  len(RefSeqMap[currentChromNO])-1+RefSeqMap[currentChromNO][0]+1 < Rearmostpos:
                    bases_to_get=Rearmostpos -len(RefSeqMap[currentChromNO])-RefSeqMap[currentChromNO][0]
                    lins_to_get=int(bases_to_get/minintervalbetweengenes_basesperfaline)+1
                    curposoffilehandler = reffa.tell()
                    reffa_suplemtry = open(options.reffa, 'r')
                    reffa_suplemtry.seek(curposoffilehandler)
                    RefSeqMap_suplemtry, currentchromNo_suplemtry,nextChromNO = Util.getRefSeqMap(reffa_suplemtry, currentChromNO=currentChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1,linesOnce=lins_to_get)
                    if currentchromNo_suplemtry!=currentChromNO:
                        print("this if block just for test")
                        exit(-1)
                    RefSeqMap[currentChromNO] += RefSeqMap_suplemtry[currentchromNo_suplemtry][1:]
                    reffa_suplemtry.close()
                    nearestGenes = Util.genes(gtfMap[currentChromNO], currentBaselocinGenome, RefSeqMap[currentChromNO])
                print(idx_RefSeq,currentBaselocinGenome,frontmostpos,Rearmostpos)
            
        else:
#            lastposofdepthfp = depthfile.depthfilefp.tell()
            if cns_string.find(">")!=-1:
                ntoseq=cns_string.find("\n",cns_string.find(">"))
                print(cns_string[0:ntoseq+1],end="",file=outcns)
                cns_string=cns_string[ntoseq+1:]
            while len(cns_string)>60:
                print(cns_string[0:60],end="\n",file=outcns)
                cns_string=cns_string[60:]
#           print(cns_string, end="", file=outcns)
            if nextChromNO == currentChromNO:
#                cns_string = ""
                RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=nextChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
            else:
                print(cns_string,end="",file=outcns)
                cns_string = "\n>" + nextChromNO + "\n"
                RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=nextChromNO, preBaseTotal=0)
        print("\t\tone loop end:", currentChromNO)
    else:
        print("finish")        
        
        
            
    
#gtffile.close()    
reffa.close()
vcffile.close()
#covfile.close()    
outcns.close()
outaa.close()
outcdscns.close()

