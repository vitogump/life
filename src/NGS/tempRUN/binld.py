'''
Created on 2015-3-7

@author: liurui
'''

from optparse import OptionParser
import re

from NGS.BasicUtil import *


parser = OptionParser()
parser.add_option("-d", "--datafilename", dest="datafilename",nargs=5,help="file_name   col_idx_to_bined1 col_idx_to_bined2 startNo endNo")
parser.add_option("-i", "--interval", dest="interval", nargs=3,help="minvalue maxvalue d_increase")

parser.add_option("-m", "--coltomean", dest="coltomean",help="coltomean")
parser.add_option("-o","--outfileprename",dest="outfileprename",help="outfilepreName with path")

                                                                                                                                                          
(options, args) = parser.parse_args()

print(options.datafilename)
print(options.interval)
minvalue=float(options.interval[0])
maxvalue = float(options.interval[1])
d_increase = float(options.interval[2])

dataFileNamePattern=options.datafilename[0]
col_to_bined1=int(options.datafilename[1])
col_to_bined2=int(options.datafilename[2])
startNo=int(options.datafilename[3])
endNo=int(options.datafilename[4])
print(startNo,endNo)
dataFileNames=tuple([re.sub(r"myNtosub",str(i),dataFileNamePattern) for i in range(startNo,endNo+1)])


col_to_mean=int(options.coltomean)
if __name__ == '__main__':
    intervalFileName=options.outfileprename+".interval"
    intervalfile=open(intervalFileName,"w")
    while minvalue + d_increase<= maxvalue :
        print(str(minvalue),str(minvalue+d_increase),sep="\t",file=intervalfile)
        minvalue+=d_increase
    else:
        if minvalue<maxvalue:
            print(str(minvalue),str(maxvalue),sep="\t",file=intervalfile)
        intervalfile.close()
    
    intervalMap_count,intervalMap_mean=Util.distributionfuncdraft(intervalFileName,dataFileNames, col_to_bined1=col_to_bined1,col_to_bined2=col_to_bined2, col_to_mean=col_to_mean)
    outfile=open(options.outfileprename,'w')
    for a,b in sorted(intervalMap_count.keys()):
        print(str(a),str(b),str(intervalMap_count[a,b]),str(intervalMap_mean[a,b]),sep="\t",file=outfile)
    outfile.close()