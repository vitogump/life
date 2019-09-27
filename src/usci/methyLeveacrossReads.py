'''
Created on 2019年9月27日

@author: liurui
'''
from optparse import OptionParser
import pysam, numpy,copy

parser = OptionParser()
parser.add_option("-b", "--bamfile", dest="bamfile",help="write report to FILE")
parser.add_option("-o","--outfileprename",dest="outfilepreName",help="outfilepreName with path")

(options, args) = parser.parse_args()

of=open(options.outfilepreName,'w')
methcounts={}
unmethcounts={}
if __name__ == '__main__':
    bamfile=pysam.AlignmentFile(options.bamfile.strip(),'rb')
    for chrom in bamfile.header['SQ']:
        for read in bamfile.fetch(chrom['SN']):
            methstats=read.get_tag('XM')
            if len(methstats) not in methcounts:
                methcounts[len(methstats)]=[0]*len(methstats)
                unmethcounts[len(methstats)]=[0]*len(methstats)
            for i in range(len(methstats)):
                if methstats[i].islower():
                    unmethcounts[len(methstats)][i]+=1
                elif methstats[i].isupper():
                    methcounts[len(methstats)][i]+=1
    for l in sorted(methcounts.keys()):
        print(*methcounts[l],sep="\t",file=of)
        print(*unmethcounts[l],sep="\t",file=of)
    bamfile.close()
            