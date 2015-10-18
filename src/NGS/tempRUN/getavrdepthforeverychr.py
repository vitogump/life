'''
Created on 2015-8-23

@author: liurui
'''
from optparse import OptionParser
import os
import re


parser = OptionParser()
parser.add_option("-b","--bamfilelist",dest="bamfilelist",action="append",help="vcftablename filerecord_allname_in_depthfiletitle_belongtothisvcfpop")
parser.add_option("-c","--chromlistfilename",dest="chromlistfilename")
parser.add_option("-o","--outfileprewithpath",dest="outfileprewithpath")
(options, args) = parser.parse_args()
if __name__ == '__main__':
    chromlistfile=open(options.chromlistfilename,"r")
    chromlist=[]
    for chrrow in chromlistfile:
        chrrowlist=re.split(r'\s+',chrrow.strip())
        chromlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))
    chromlistfile.close()
    for bamfilelistfilenamewithpath in options.bamfilelist:
        bamfilelistfile=open(bamfilelistfilenamewithpath,"r")
        for bamfilename in bamfilelistfile:
            bamname=re.search(r"[^/]+$",bamfilename.strip()).group(0)
            chromlistfilenamewithoutpath=re.search(r"[^/]+$",options.chromlistfilename).group(0)
            outfile=open(options.outfileprewithpath+chromlistfilenamewithoutpath+bamname+"_depForChr","w")
            for currentID,currentchrLen in chromlist:
                if not os.path.isfile(bamfilename.strip()+".bai"):
                    nnnnnnn=re.search(r"(.*)[.]bam",bamfilename.strip()).group(1)
                    print("cp "+nnnnnnn+".bai "+bamfilename.strip()+".bai")
                    os.system("cp "+nnnnnnn+".bai "+bamfilename.strip()+".bai")
                print("/pub/tool/samtools-0.1.19/bin/samtools depth -r "+currentID+" "+bamfilename.strip()+" | awk '{sum += $3} END {print sum / NR}'")
                f=os.popen("/pub/tool/samtools-0.1.19/bin/samtools depth -r "+currentID+" "+bamfilename.strip()+" | awk '{sum += $3} END {print sum / NR}'")
                meandepth=f.readline().strip()
                print(currentID,currentchrLen,meandepth,sep="\t",file=outfile)
                f.close()
            outfile.close()
        bamfilelistfile.close()