'''
Created on 2019年8月12日

@author: RuiLiu
'''
from optparse import OptionParser

import pysam, numpy

parser = OptionParser()
parser.add_option("-b", "--bamfile", dest="bamfile",help="write report to FILE")
parser.add_option("-r", "--reffa", dest="reffa",help="far sure but with few locs")
parser.add_option("-o","--outfileprename",dest="outfilepreName",help="outfilepreName with path")
parser.add_option("-i", "--interval", dest="interval", nargs=3, help="minvalue maxvalue breaks. divid the delta (P1Freq-P2Freq)")
parser.add_option("-m","--mindepth",dest="mindepth",help="mindepth for both archicpop and ancestralallel")#
# parser.add_option("-e", "--bedfile", action="append",dest="bedfile",default=[],help="measure region of bedfile")
parser.add_option("-s","--speciesesName",action="append",dest="speciesesName",default=[],help="speciesName in table")#
                                                                                                                                                          
(options, args) = parser.parse_args()
chromMap={
    "chrM":16571,
"chr1":249250621,
"chr2":243199373,
"chr3":198022430,
"chr4":191154276,
"chr5":180915260,
"chr6":171115067,
"chr7":159138663,
"chr8":146364022,
"chr9":141213431,
"chr10":135534747,
"chr11":135006516,
"chr12":133851895,
"chr13":115169878,
"chr14":107349540,
"chr15":102531392,
"chr16":90354753,
"chr17":81195210,
"chr18":78077248,
"chr19":59128983,
"chr20":63025520,
"chr21":48129895,
"chr22":51304566,
"chrX":155270560,
"chrY":59373566}

minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
dincrease = (maxvalue - minvalue) / breaks
delta_DerAftotal={}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    delta_DerAftotal[(minvalue,minvalue + dincrease)]=0
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
def Compute_GC(seq):
    GC=0
    total=0
    for i in range(0,len(seq)):
        if seq[i].upper() in ["A","T","G","C"]:
            total=total+1
        if seq[i].upper() in ["G","C"]:
            GC=GC+1
    return float(GC)/total*100
if __name__ == '__main__':
    bamfile=pysam.Samfile(options.bamfile.strip(),'rb')
    print(bamfile.references)
#     bamfile=pysam.AlignmentFile(options.bamfile.strip(),'rb')
    for chrom in bamfile.references:
        print(chrom)
        for read  in bamfile.fetch(chrom):
            rseq=read.query_sequence
            GCcontent=Compute_GC(rseq)
            for a,b in sorted(delta_DerAftotal.keys()):
                if GCcontent >a and GCcontent<=b:
                    delta_DerAftotal[(a,b)]+=1
                    break
        print(delta_DerAftotal)
        
    for a,b in sorted(delta_DerAftotal.keys()):
        print(a,b,delta_DerAftotal[a,b],file=)
    bamfile.close()
        