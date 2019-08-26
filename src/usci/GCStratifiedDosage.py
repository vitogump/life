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
Zoutfile=open(options.outfilepreName+"_normlized",'w')
scaledoutfile=open(options.outfilepreName+"_scaled",'w')
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
bins_arr=[[],[]]#lower than minvalue, more than maxvalue
exception=[[],[]];std1=[[],[]]
while minvalue<=maxvalue - dincrease :
    print(minvalue,minvalue + dincrease)
    bin_GCstratified_template[(minvalue,minvalue + dincrease)]=0
    bin_GCstratified_mean[(minvalue,minvalue + dincrease)]=0
    bins_arr.append([]);exception.append([]);std1.append([])# each elem is a list of genome bin
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
#         if chrom['SN']=="chr2":
#             break
        for str_bp in range(0,chrom['LN'],100000):# for every genomebin
            bin_GCstratified_colect=copy.deepcopy(bin_GCstratified_template)
            for GCbin_index in range(len(bins_arr)):
                bins_arr[GCbin_index].append(0) 
            print(str_bp,str_bp+100000)
            for read  in bamfile.fetch(chrom['SN'],str_bp,str_bp+100000):
                rseq=read.query_sequence
                GCcontent=Compute_GC(rseq)
                GCbin_index=0
                for a,b in sorted(bin_GCstratified_colect.keys()):
                    if GCcontent >a and GCcontent<=b:
                        bin_GCstratified_colect[(a,b)]+=1
                        bins_arr[GCbin_index][-1]=bin_GCstratified_colect[(a,b)]
                        bin_GCstratified_mean[(a,b)]+=1
                        GCbin_index+=1
                        break
                    elif GCcontent <=a:
                        bins_arr[GCbin_index][-1]+=1
                    elif GCcontent >b:
                        bins_arr[GCbin_index][-1]+=1
                    GCbin_index+=1
            bins_GCstratified.append(bin_GCstratified_colect)
    print(len(bins_GCstratified))
    
#     for genomeBin in bins_GCstratified:
#         for a,b in sorted(genomeBin.keys()):
    for GCbin_index in range(len(bins_arr)):
        exception[GCbin_index]=numpy.mean(bins_arr[GCbin_index])
        std1[GCbin_index]=numpy.std(bins_arr[GCbin_index],ddof=1)
    print(exception,std1)
    for a,b in sorted(bin_GCstratified_mean.keys()):
        bin_GCstratified_mean[a,b]=bin_GCstratified_sum[a,b]/len(bins_GCstratified)#overwrite the bin_GCstratified_mean which is actually a overlapped sum
    #bin_GCstratified_mean ≈ exception
    print("binstr\tbinend","<",end="\t",file=Zoutfile);print("binstr\tbinend","<",end="\t",file=scaledoutfile)
    for a,b in sorted(bin_GCstratified_colect.keys()):
        print(str(a)+"-"+str(b),end="\t",file=Zoutfile)
        print(str(a)+"-"+str(b),end="\t",file=scaledoutfile)
    print(">",file=Zoutfile);print(">",file=scaledoutfile)
    idx=0
    for bin_GCstratified_colect in bins_GCstratified:
        print(idx*100000,(idx+1)*100000,sep="\t",end="\t",file=Zoutfile)
        print(idx*100000,(idx+1)*100000,sep="\t",end="\t",file=scaledoutfile)
        for GCidx in range(len(bins_arr)):
            
            try:
                print(bins_arr[GCidx][idx]/exception[GCidx],end="\t",file=scaledoutfile)# use bin_GCstratified_mean maybe better
                print((bins_arr[GCidx][idx]-exception[GCidx])/std1[GCidx],end="\t",file=Zoutfile)
            except:
                print("NA",end="\t",file=scaledoutfile)
                print("NA",end="\t",file=Zoutfile)
        print("",file=scaledoutfile);print("",file=Zoutfile)
        idx+=1
#         for a,b in sorted(bin_GCstratified_colect.keys()):    
#             print(a,b,bin_GCstratified_sum[a,b],sep="\t",file=Zoutfile)
#             print(a,b,bin_GCstratified_sum[a,b],sep="\t",file=scaledoutfile)
    bamfile.close();scaledoutfile.close();Zoutfile.close()
        