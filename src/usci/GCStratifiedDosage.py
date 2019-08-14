'''
Created on 2019年8月12日

@author: RuiLiu
'''
from optparse import OptionParser

import pysam, numpy,copy

parser = OptionParser()
parser.add_option("-b", "--bamfile", dest="bamfile",help="write report to FILE")
parser.add_option("-r", "--reffa", dest="reffa",help="far sure but with few locs")
parser.add_option("-o","--outfileprename",dest="outfilepreName",help="outfilepreName with path")
parser.add_option("-i", "--interval", dest="interval", nargs=3, help="minvalue maxvalue breaks. divid the delta (P1Freq-P2Freq)")
parser.add_option("-m","--mindepth",dest="mindepth",help="mindepth for both archicpop and ancestralallel")#
# parser.add_option("-e", "--bedfile", action="append",dest="bedfile",default=[],help="measure region of bedfile")
parser.add_option("-s","--speciesesName",action="append",dest="speciesesName",default=[],help="speciesName in table")#
                                                                                                                                                          
(options, args) = parser.parse_args()
outfile=open(options.outfilepreName,'w')
minvalue = float(options.interval[0])
maxvalue = float(options.interval[1])
breaks = int(options.interval[2])
dincrease = (maxvalue - minvalue) / breaks
bin_GCstratified_template={}#use this as template 
bin_GCstratified_mean={}# recod sum value
bin_GCstratified_sum={
(30.0,34.0): 4249574,
(34.0,38.0): 4917351,
(38.0,42.0): 4738324,
(42.0,46.0): 3974121,
(46.0,50.0): 3214421,
(50.0,54.0): 2718319,
(54.0,58.0): 2300071,
(58.0,62.0): 1574971,
(62.0,66.0): 802171,
(66.0,70.0): 336057}
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    bin_GCstratified_template[(minvalue,minvalue + dincrease)]=0
    bin_GCstratified_mean[(minvalue,minvalue + dincrease)]=0
#     delta_DerAf[(minvalue,minvalue + dincrease)]={"BinP1andP2":[],"BinP1orP2":[],"BBAA":[],"BABA":[],"ABBA":[]}
    minvalue += dincrease
bins_GCstratified=[]
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
#     bamfile=pysam.Samfile(options.bamfile.strip(),'rb')
    bamfile=pysam.AlignmentFile(options.bamfile.strip(),'rb')
    print(len(bamfile.header['SQ']),len(bamfile.references))
    
    for chrom in bamfile.header['SQ']:
        print(chrom,chrom['SN'])
        for str_bp in range(0,chrom['LN'],100000):
            print(str_bp,str_bp+100000)
            bin_GCstratified_colect=copy.deepcopy(bin_GCstratified_template)
            for read  in bamfile.fetch(chrom['SN'],str_bp,str_bp+100000):
                rseq=read.query_sequence
                GCcontent=Compute_GC(rseq)
                for a,b in sorted(bin_GCstratified_colect.keys()):
                    if GCcontent >a and GCcontent<=b:
                        bin_GCstratified_colect[(a,b)]+=1
                        bin_GCstratified_mean[(a,b)]+=1
                        break
            bins_GCstratified.append(bin_GCstratified_colect)
    print(len(bins_GCstratified))
    
    
    for a,b in sorted(bin_GCstratified_mean.keys()):
        bin_GCstratified_mean[a,b]=bin_GCstratified_sum[a,b]/len(bins_GCstratified)#overwrite the bin_GCstratified_mean which is actrually a overlaped sum
    for bin_GCstratified_colect in bins_GCstratified:
        for a,b in sorted(bin_GCstratified_colect.keys()):    
            print(a,b,bin_GCstratified_sum[a,b],sep="\t",file=outfile)
    bamfile.close();outfile.close()
        