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
if len(sys.argv) != 6:
    print(len(sys.argv))
    print("python GetConveredConsensusSeq.py [ref.fa] [gtffile] [vcffile(withheader)] [genomeCoveragefile] [outfileprename] ")
    exit(-1)
parser = OptionParser()
parser.add_option("-r", "--reffa", dest="reffa",
                  help="reference.fa", metavar="FILE")
parser.add_option("-g", "--gtffile", dest="gtffile", help="gtffile")
parser.add_option("-v", "--vcffile", dest="variant", help="variants")
parser.add_option("-d", "--genomeCoveragefile", dest="genomedepth", help="genomeCoveragefile")
parser.add_option("-o", "--outfileprename", dest="outfileprename", help="outfileprename")
parser.add_option("-C", "--whetherwriteconsensus", dest="cnsornot", action='store_false', help="use this option means you will not print the cns.fa file")
parser.add_option("-n", "--speciesname", dest="species", help="species name")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
                                                                                                                                                          
(options, args) = parser.parse_args()
reffa = open(options.reffa, 'r')
gtffile = open(options.gtffile, 'r')
vcffile = open(options.variants, 'r')
covfile = open(options.genomedepth, 'r')
cns_string=""
aa_string=""
cdscns_string=""
outcns = open(options.outfileprename + "_cns.fa", 'w')
outaa = open(options.outfileprename + "_aa.fa", 'w')
outcdscns = open(options.outfileprename + "_cdscns.fa", 'w')
if __name__ == '__main__':
    depthfile = Util.GATK_depthfile(options.genomedepth, options.genomedepth + ".index")
    species_idx = depthfile.title.index("Depth_for_" + options.species)
    vcfpop = VCFutil.VCF_Data(options.variants)  # new a class
    RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa)
#    Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=currentchr, preBaseTotal=RefSeqMap[lastchromNo][0] + len(RefSeqMap[lastchromNo]) - 1, linesOnce)
    lastposofdepthfile = 0#because this time RefSeqMap[0] is 0
    while nextChromNO != "end of the reffile":
        currentBaselocinGenome = RefSeqMap[0]
        statue = depthfile.set_depthfilehandler(currentChromNO, currentBaselocinGenome, lastposofdepthfile)
        for currentBase in RefSeqMap[currentChromNO][1:]:#base by base
            currentBaselocinGenome += 1
            depth_chrom, depth_pos, depth_line, depth_linelist = depthfile.getnextposline()
            if depth_chrom != currentChromNO or depth_pos!=currentBaselocinGenome:
                statue = depthfile.set_depthfilehandler(currentChromNO, currentBaselocinGenome, lastposofdepthfile)
                if statue =="didn't find":
                    print("warning:"+options.genomedepth+"didn't have this genome postion:"+currentChromNO+" , "+currentBaselocinGenome)
                    cns_string+="n"
            else:
                if int(depth_linelist[species_idx]) >= mindepth:
                    pass
                else:
                    cns_string+="n"
        lastposofdepthfile = depthfile.depthfilehandler.tell()
        print(cns_string,end="",file=outcns)
        if nextChromNO == currentChromNO:
            cns_string=""
        else:
            cns_string="/n>"+nextChromNO+"/n"
        RefSeqMap, currentChromNO, nextChromNO = Util.getRefSeqMap(refFastafilehander=reffa, currentChromNO=nextChromNO, preBaseTotal=RefSeqMap[currentChromNO][0] + len(RefSeqMap[currentChromNO]) - 1)
            
    
    
    
